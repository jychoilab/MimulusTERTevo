#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


# --------------------------------------------------
# INPUT
# --------------------------------------------------

infile = "TRAP_read_counts.tsv"


# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------

df = pd.read_csv(
    infile,
    sep="\t"
)


# --------------------------------------------------
# Normalize TRAP products
# --------------------------------------------------

df["TR1_plus_TR2"] = (
    (df["TR1_reads"] + df["TR2_reads"]) /
    df["Combined_E_coli_coverage"]
)


df["TR1"] = (
    df["TR1_reads"] /
    df["Combined_E_coli_coverage"]
)


df["TR2"] = (
    df["TR2_reads"] /
    df["Combined_E_coli_coverage"]
)


df.to_csv(
    "TRAP_normalized_input.tsv",
    sep="\t",
    index=False
)



# --------------------------------------------------
# Species
# --------------------------------------------------

df["Species"] = (
    df["Sample"]
    .str.extract("(car|lew|ver)")
)


species_order = [
    "car",
    "lew",
    "ver"
]


species_labels = [
    "Cardinalis",
    "Lewisii",
    "Verbenaceus"
]



# --------------------------------------------------
# Statistics
# --------------------------------------------------

def p_to_star(p):

    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"



def add_sig(ax, x1, x2, y, text, color="black"):

    h = 0.08

    ax.plot(
        [x1, x1, x2, x2],
        [y, y+h, y+h, y],
        linewidth=1.5,
        color=color
    )


    ax.text(
        (x1+x2)/2,
        y+h,
        text,
        ha="center",
        va="bottom",
        fontsize=11,
        color=color
    )



# --------------------------------------------------
# Plot 1: TR1 + TR2 combined
# --------------------------------------------------

def plot_TR1_plus_TR2(with_stars):


    temp = df.copy()

    temp["log10"] = np.log10(
        temp["TR1_plus_TR2"]
    )


    means = []
    sds = []


    for sp in species_order:

        values = temp[
            temp["Species"] == sp
        ]["log10"]

        means.append(values.mean())
        sds.append(values.std())



    fig, ax = plt.subplots(
        figsize=(7,6)
    )


    x = np.arange(3)


    ax.bar(
        x,
        means,
        yerr=sds,
        capsize=5,
        color="lightgray",
        edgecolor="black",
        linewidth=1.8
    )


    # replicate points

    for i,sp in enumerate(species_order):

        values = temp[
            temp["Species"] == sp
        ]["log10"]


        ax.scatter(
            np.repeat(i,len(values)),
            values,
            color="black",
            s=35,
            zorder=3
        )


    ax.set_xticks(x)

    ax.set_xticklabels(
        species_labels
    )


    ax.set_xlabel("")


    ax.set_ylabel(
        "Normalized TRAP products\n(log10 reads per 1× bacterial genome coverage)",
        fontsize=12
    )


    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


    if with_stars:


        ymax=max(means)+max(sds)


        for k in range(3):

            for j in range(k+1,3):

                g1=temp[
                    temp.Species==species_order[k]
                ]["log10"]


                g2=temp[
                    temp.Species==species_order[j]
                ]["log10"]


                _,p=mannwhitneyu(
                    g1,
                    g2,
                    alternative="two-sided"
                )


                add_sig(
                    ax,
                    k,
                    j,
                    ymax+k*0.25,
                    p_to_star(p)
                )


    plt.tight_layout()


    outfile = (
        "TR1_plus_TR2_bar_log10_stars.pdf"
        if with_stars
        else
        "TR1_plus_TR2_bar_log10_clean.pdf"
    )


    plt.savefig(
        outfile,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()



# --------------------------------------------------
# Plot 2: TR1 vs TR2 grouped bars
# --------------------------------------------------

def plot_TR1_vs_TR2(with_stars):


    long=df.melt(

        id_vars=[
            "Sample",
            "Species"
        ],

        value_vars=[
            "TR1",
            "TR2"
        ],

        var_name="Condition",

        value_name="value"
    )


    long["log10"]=np.log10(
        long["value"]
    )


    colors={
        "TR1":"#ADD8E6",
        "TR2":"#800080"
    }


    fig,ax=plt.subplots(
        figsize=(8,6)
    )


    width=0.35


    x=np.arange(3)



    for idx,condition in enumerate(["TR1","TR2"]):


        means=[]
        sds=[]


        for sp in species_order:


            values=long[
                (long.Species==sp) &
                (long.Condition==condition)
            ]["log10"]


            means.append(values.mean())

            sds.append(values.std())



        shift=(-width/2 if condition=="TR1" else width/2)



        ax.bar(

            x+shift,

            means,

            width,

            yerr=sds,

            capsize=4,

            color=colors[condition],

            edgecolor="black",

            linewidth=1.5,

            label=condition
        )


        # points

        for i,sp in enumerate(species_order):


            values=long[
                (long.Species==sp) &
                (long.Condition==condition)
            ]["log10"]


            ax.scatter(

                np.repeat(
                    x[i]+shift,
                    len(values)
                ),

                values,

                color="black",

                s=30,

                zorder=3
            )



    ax.set_xticks(x)

    ax.set_xticklabels(
        species_labels
    )


    ax.set_xlabel("")


    ax.set_ylabel(
        "Normalized TRAP products\n(log10 reads per 1× bacterial genome coverage)",
        fontsize=12
    )


    ax.legend(
        frameon=False
    )


    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)



    if with_stars:


        ymax=long["log10"].max()+0.3


        for condition,color in colors.items():


            shift=(
                -width/2
                if condition=="TR1"
                else
                width/2
            )


            for i in range(3):

                for j in range(i+1,3):


                    g1=long[
                        (long.Species==species_order[i]) &
                        (long.Condition==condition)
                    ]["log10"]


                    g2=long[
                        (long.Species==species_order[j]) &
                        (long.Condition==condition)
                    ]["log10"]


                    _,p=mannwhitneyu(
                        g1,
                        g2,
                        alternative="two-sided"
                    )


                    add_sig(
                        ax,
                        i+shift,
                        j+shift,
                        ymax,
                        p_to_star(p),
                        color
                    )

                    ymax+=0.25



    plt.tight_layout()


    outfile = (
        "TR1_vs_TR2_bar_log10_stars.pdf"
        if with_stars
        else
        "TR1_vs_TR2_bar_log10_clean.pdf"
    )


    plt.savefig(
        outfile,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()



# --------------------------------------------------
# RUN
# --------------------------------------------------

plot_TR1_plus_TR2(False)
plot_TR1_plus_TR2(True)

plot_TR1_vs_TR2(False)
plot_TR1_vs_TR2(True)


print("Finished generating TRAP bar plots.")
