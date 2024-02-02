# open xlsx file
import pandas as pd

import csv
import os
import mygene

# open xlsx file media-1.xlsx
df = pd.read_excel("media-1.xlsx", sheet_name="S4", header=[0, 1])
print(df.columns)


# MultiIndex([('Level-1 annotation',      'p_val'),
#             ('Level-1 annotation', 'avg_log2FC'),
#             ('Level-1 annotation',      'pct.1'),
#             ('Level-1 annotation',      'pct.2'),

# #get level-1 annotation
# load renaming from csv
renaming = {}
if os.path.exists("renaming.csv"):
    with open("renaming.csv", "r") as f:
        reader = csv.reader(f)
        renaming = dict(reader)
print("using renaming", renaming)

# add manual RP11-620J15.3 to ENSG00000257698
renaming[
    "RP11-620J15.3"
] = "GIHCG"  # https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000257698;r=12:58325311-58329958
renaming["RP11-1143G9.4"] = "AC020656.1" #https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000257764;r=12:69747273-69748005;t=ENST00000548900



# renaming = {"MIR9-1HG": "C1orf61"}

unknown = []

for level_type in [
    "Level-1 annotation",
    "Level-2 annotation",
    "Level-3 annotation",
    "Level-4 annotation",
]:
    level = df[level_type].copy()
    print(level.columns)
    # remove nans
    level = level.dropna()

    # get clusters from level-1 annotation
    clusters = level["cluster"]
    genes = level["gene"]
    # print(clusters)
    # print(genes)

    # open  /Users/carbs/Downloads/Spatial_Meta/visium/A-AK40374-AK40375_225GY5LT3/outs/filtered_feature_bc_matrix/features.tsv
    # get gene names
    # get ensembl ids

    # mapp = pd.read_csv(
    #     "/Users/carbs/Downloads/Spatial_Meta/visium/A-AK40374-AK40375_225GY5LT3/outs/filtered_feature_bc_matrix/features.tsv",
    #     sep="\t",
    #     header=None,
    # )

    # /Users/carbs/Downloads/Spatial_Meta/visium/A-AK40374-AK40375_225GY5LT3/SPATIAL_RNA_COUNTER_CS/SPATIAL_RNA_COUNTER/PARSE_TARGET_FEATURES/fork0/chnk0-u63f66a0ec3/files/probe_set.csv
    # it has comments on fist 4 lines so skip them
    # mapp2= pd.read_csv("/Users/carbs/Downloads/Spatial_Meta/visium/A-AK40374-AK40375_225GY5LT3/SPATIAL_RNA_COUNTER_CS/SPATIAL_RNA_COUNTER/PARSE_TARGET_FEATURES/fork0/chnk0-u63f66a0ec3/files/probe_set.csv", skiprows=5)
    # /Users/carbs/Downloads/Spatial_Meta/visium/A-AK40374-AK40375_225GY5LT3/SPATIAL_RNA_COUNTER_CS/SPATIAL_RNA_COUNTER/_TARGETED_ANALYZER/CALCULATE_TARGETED_METRICS/fork0/join-u63f66a1acf/files/per_feature_metrics_csv.csv
    map2 = pd.read_csv(
        "/Users/carbs/Downloads/Spatial_Meta/visium/A-AK40374-AK40375_225GY5LT3/SPATIAL_RNA_COUNTER_CS/SPATIAL_RNA_COUNTER/_TARGETED_ANALYZER/CALCULATE_TARGETED_METRICS/fork0/join-u63f66a1acf/files/per_feature_metrics_csv.csv"
    )
    # extract only column with name probe_id
    # mapp2 = mapp2[["probe_id"]]

    # #split probe_ids ENSG00000000003|TSPAN6|41ef80c
    # mapp2 = mapp2["probe_id"].str.split("|", expand=True)
    # #rename columns
    # mapp2.columns = ["gene_id", "gene_name", "somehashoridk"]

    # keep only colums with name feature_type,feature_id,feature_name,
    mapp2 = map2[["feature_type", "feature_id", "feature_name"]]
    # rename columns
    mapp2.columns = ["feature_type", "gene_id", "gene_name"]

    # transfrom  b'ENSG00000130203' to ENSG00000130203 by removing b' and '
    mapp2["gene_id"] = mapp2["gene_id"].str[2:-1]

    # print(mapp2.head())

    # find TUBA1A CCL3L3

    # print(mapp2[mapp2["gene_name"] == "TUBA1A"].empty)

    # name columns
    # mapp.columns = ["ensembl", "gene", "type"]
    # print(mapp.head())

    # not found
    not_found = []

    # manal rename MIR9-1HG to C1orf61

    # map
    for gene in genes:
        # print(gene)
        # print(mapp[mapp["gene"] == gene])

        # check if gene is in renaming
        if gene in renaming:
            gene = renaming[gene]
            print("renamed", gene)

        if mapp2[mapp2["gene_name"] == gene].empty:
            print("gene not found", gene)
            not_found.append(gene)

    mg = mygene.MyGeneInfo()

    # serch missing
    query_results = mg.querymany(
        not_found, scopes="symbol", fields="ensembl.gene", species="human"
    )
    # Process query results
    for result in query_results:
        # print(result)
        # get ens id
        ensembl_ids = result.get("ensembl", {})
        print(result)
        if isinstance(ensembl_ids, list):
            ensembl_id = ensembl_ids[0]["gene"]
            print("multiple ensembl ids found for", result.get("query"), ensembl_ids)
            for ensembl_id in ensembl_ids:
                # print(ensembl_id["gene"])
                found = False
                if ensembl_id["gene"] in mapp2["gene_id"].values:
                    print("found in mapp", ensembl_id["gene"])
                    # add to renaming
                    renaming[result.get("query")] = mapp2[
                        mapp2["gene_id"] == ensembl_id["gene"]
                    ]["gene_name"].tolist()[0]
                    print("renaming changed", renaming)
                    found = True
                    break
            if not found:
                print("not found in mapp", result.get("query"), ensembl_ids)
                unknown.append(result.get("query"))

        # cehck
        elif isinstance(ensembl_ids, dict):
            # {'query': 'RP11-620J15.3', 'notfound': True}
            # check
            if result.get("notfound"):
                print("not found in mapp", result.get("query"))
                unknown.append(result.get("query"))
                continue

            ensembl_id = ensembl_ids["gene"]
            # check if ensembl id is in mapp
            if ensembl_id in mapp2["gene_id"].values:
                print("found in mapp", ensembl_id)
                # add to renaming
                renaming[result.get("query")] = mapp2[mapp2["gene_id"] == ensembl_id][
                    "gene_name"
                ].tolist()[0]
                print("renaming changed", renaming)
            else:
                print("not found in mapp", result.get("query"), ensembl_id)
                unknown.append(result.get("query"))

        else:
            print("nohit", result.get("query"))
            unknown.append(result.get("query"))

    # Create a MyGeneInfo object

    # map gene to gene name and add gene_ens

    # Assuming 'level_1' is a DataFrame and 'gene' is the column with gene symbols
    gene_symbols = level[
        "gene"
    ].unique()  # Get unique gene symbols to avoid duplicate queries

    # group by cluster and get genes
    grouped = level.groupby("cluster")
    # print(grouped.groups)
    print("------------------")

    # open csv it has colums named List,Name,ID,Feature Type

    for cluster, group in grouped:
        # print(cluster)

        # find gene_id for each gene name
        # get gene name
        gene_names = group["gene"].tolist()
        # get gene id check if it is in mapp2 else skip
        gene_ids = []
        for gene_name in gene_names:
            # check if gene is in renaming
            if gene_name in renaming:
                gene_name = renaming[gene_name]
                print("renamed", gene_name)

            gene_id = mapp2[mapp2["gene_name"] == gene_name]["gene_id"].tolist()
            if gene_id:
                gene_ids.append(gene_id[0])
            else:
                print(
                    "gene not found",
                    gene_name,
                    "in cluster",
                    cluster,
                    " Dropping this feature",
                )
                gene_ids.append("not found")
                unknown.append(gene_name)

        # add gene_id to group
        group["gene_id"] = gene_ids

        # add these to level_1
        level.loc[group.index, "gene_id"] = gene_ids
        # set not found to nan
        level.loc[level["gene_id"] == "not found", "gene_id"] = None

        # drop rows with not found
        group = group[group["gene_id"] != "not found"]

        # print(group)

    # print(level)

    # get NANs
    nan_rows = level[level["gene_id"].isnull()]
    # print(nan_rows)

    # drop nan rows
    level = level.dropna()
    # print(level)

    # aduplicate cluster column and rename ti to List
    level["List"] = level["cluster"]
    # duplicate gene column and rename it to Name
    level["Name"] = level["gene"]
    # duplicate gene_id column and rename it to ID
    level["ID"] = level["gene_id"]

    # add FeatureType column and have all values as Gene Expression
    level["Feature Type"] = "Gene Expression"

    # save to csv only List,Name,ID,Feature Type
    level[["List", "Name", "ID", "Feature Type"]].to_csv(
        level_type + ".csv", index=False
    )

# save renaming
with open("renaming.csv", "w") as f:
    w = csv.writer(f)
    w.writerows(renaming.items())

# save unknown
with open("unknown_genes.txt", "w") as f:
    for item in unknown:
        f.write("%s\n" % item)
