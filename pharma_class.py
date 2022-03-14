import pandas as pd
import numpy as np
import sys
import os

# paths
reference_path = sys.path[0] + "/Configuration"
recommendations_path = reference_path + "/Dosage_Recommendations_December_2021.csv"
alternative_recommendations_path = reference_path + "/Alternative_Recommendations.xls"
snp_path = reference_path + "/snp_final.txt"
star_path = reference_path +"/star_final.txt"

# diplotype phenotypes
diplotype_phenotype = reference_path + "/Diplotype_Phenotype_December_2021/All_December_2021.csv"

#exceptions
exceptions_other_names_path = reference_path + "/Diplotype_Phenotype_December_2021/Exceptions_Other_Names.csv"
exceptions_rsid_path = reference_path + "/Diplotype_Phenotype_December_2021/Exceptions_RSID_December_2021.csv"

# load tables
recommendations = pd.read_csv(recommendations_path, dtype="object")
alternative_recommendations = pd.read_excel(alternative_recommendations_path, dtype="object")
snp = pd.read_csv(snp_path, sep="\t", dtype="object")
star = pd.read_csv(star_path, sep="\t", dtype="object")
all = pd.read_csv(diplotype_phenotype, dtype="object")
exceptions_other_names = pd.read_csv(exceptions_other_names_path,dtype="object")
exceptions_rsid = pd.read_csv(exceptions_rsid_path,dtype="object")

# gene lists

# special considerations for phenotype type
exceptions_gene_other_names = ["RYR1", "CACNA1S", "CFTR", "G6PD", "CYP2C9", "DPYD"]
exceptions_gene_rsid = ["IFNL3", "VKORC1", "CYP4F2"]

# other_names
other_names = ["RYR1", "CACNA1S", "CFTR", "G6PD", "IFNL3", "VKORC1", "DPYD"]

# warfarin
warfarin_genes = ["CYP2C9","CYP4F2", "VKORC1"]

pd.set_option('display.max_columns', None)

class SubjectPharmaBuilder:
    def __init__(self, stargazer, rsid):

        self.raw = pd.read_csv(stargazer, sep="\t", dtype="object")
        filesize = os.path.getsize(rsid)
        if filesize != 0:
            self.rsid = pd.read_csv(rsid, sep="\t", dtype="object", header=None)
        else:
            self.rsid = pd.DataFrame()

    def __add_other_names(self):
        raw = self.raw
        other_names = star[["gene", "name", "other_names"]]
        raw["GENE"] = raw["gene"].apply(lambda x: x.upper())

        # prep data for merging
        other_names = other_names.rename(columns={"name":"hap1_main"})
        raw = raw.merge(other_names, on=["gene", "hap1_main"], how="left")
        raw = raw.rename(columns={"other_names":"other_names1"})

        # prep data for merging
        other_names = other_names.rename(columns={"hap1_main":"hap2_main"})
        raw = raw.merge(other_names, on=["gene", "hap2_main"], how="left")
        raw = raw.rename(columns={"other_names":"other_names2"})

        # save other_names into output
        self.raw = raw

        return self

    def __add_rsid(self):
        raw = self.raw
        rsid = snp[["wgene", "pos", "wt", "var", "id"]]

        # prep data for merging
        rsid = rsid.rename(columns={"wgene":"gene", "pos":"pos1","wt": "wt1", "var":"var1", "id":"rsid1"})
        raw.hap1_main_core = raw.hap1_main_core.apply(lambda x: x[1:len(x) - 1].split(":"))
        raw["pos1"] = [x[0] if len(x) > 1 else "." for x in raw.hap1_main_core]
        raw["wt1"] = [x[1].split(">")[0] if len(x)> 1 else "." for x in raw.hap1_main_core]
        raw["var1"] = [x[1].split(">")[1] if len(x)> 1 else "." for x in raw.hap1_main_core]
        raw = raw.merge(rsid, how="left").fillna("")

        # prep data for merging
        rsid = rsid.rename(columns={"pos1":"pos2", "wt1":"wt2", "var1":"var2", "rsid1":"rsid2"})
        raw.hap2_main_core = raw.hap2_main_core.apply(lambda x: x[1:len(x) - 1].split(":"))
        raw["pos2"] = [x[0] if len(x) > 1 else "." for x in raw.hap2_main_core]
        raw["wt2"] = [x[1].split(">")[0] if len(x)> 1 else "." for x in raw.hap2_main_core]
        raw["var2"] = [x[1].split(">")[1] if len(x)> 1 else "." for x in raw.hap2_main_core]
        raw = raw.merge(rsid, how="left").fillna("")

        self.raw = raw
        return self


    def __get_pharmgkb_phenotypes(self):
        raw = self.raw

        # normal cases except DPYD
        other = all
        other = other.rename(columns={"Variant 1":"hap1_main", "Variant 2": "hap2_main"})
        raw = raw.merge(other, how="left", on=["GENE", "hap1_main", "hap2_main"])

        # exception cases
        for gene in exceptions_gene_other_names:

            # extract gene index from raw
            gene_index = raw.loc[raw.GENE == gene].index[0]
            variant1, variant2 = raw.loc[raw.GENE == gene]["hap1_main"].item(), raw.loc[raw.GENE == gene]["hap2_main"].item()
            other_names1, other_names2 = raw.loc[raw.GENE == gene]["other_names1"].item(), raw.loc[raw.GENE == gene]["other_names2"].item()
            if gene != "DPYD":
                variant_list = exceptions_other_names.loc[exceptions_other_names.GENE == gene]["Variant"].values

            #RYR1 and CACNA1S
            if gene == "RYR1" or gene == "CACNA1S":
                if other_names1 in variant_list or other_names2 in variant_list:
                    raw.at[gene_index, "Phenotype"] = "Malignant Hyperthermia Susceptible"
                else:
                    raw.at[gene_index, "Phenotype"] = "Uncertain Susceptibility"

            #CFTR
            if gene == "CFTR":
                if other_names1 in variant_list or other_names2 in variant_list:
                    if "F508del" in [variant1,other_names1] and "F508del" in [variant2,other_names2]:
                        raw.at[gene_index, "Phenotype"] = "Ivacaftor non-responsive in CF patients"
                    else:
                        raw.at[gene_index, "Phenotype"] = "Ivacaftor responsive in CF patients"
                else:
                    raw.at[gene_index, "Phenotype"] = "Ivacaftor non-responsive in CF patients"

            if gene == "G6PD":
                variant_class = exceptions_other_names.loc[exceptions_other_names.GENE == gene]

                class_iv_variants = variant_class.loc[variant_class.Class == "IV"]["Variant"].values
                num_class_iv = 0
                num_class_other = 0

                for variant in [other_names1,other_names2]:
                    if variant in class_iv_variants:
                        num_class_iv += 1
                    if variant in variant_list:
                        num_class_other += 1

                if num_class_iv > 0:
                    raw.at[gene_index, "Phenotype Male"] = "Normal Metabolizer"
                    if num_class_iv == 2:
                        raw.at[gene_index, "Phenotype Female"] = "Normal Metabolizer"
                    else:
                        raw.at[gene_index, "Phenotype Female"] = "Variable"
                elif num_class_other > 0:
                    raw.at[gene_index, "Phenotype Male"] = "Deficient"
                    if num_class_other == 2:
                        raw.at[gene_index, "Phenotype Female"] = "Deficient"
                    else:
                        raw.at[gene_index, "Phenotype Female"] = "Variable"
                else:
                    raw.at[gene_index, "Phenotype Male"], raw.at[gene_index, "Phenotype Female"] = "Unknown",  "Unknown"

                raw.at[gene_index, "Phenotype"] = "Female:\n" + raw.at[gene_index, "Phenotype Female"] + "\nMale:\n" + raw.at[gene_index, "Phenotype Male"]

            if gene == "CYP2C9":
                raw.at[gene_index, "Phenotype Warfarin"] = "Normal Metabolizer"
                if variant1 in variant_list:
                    if variant1 in ["*2", "*3"]:
                        raw.at[gene_index, "Phenotype Warfarin"] = "Poor Function"
                    else:
                        raw.at[gene_index, "Phenotype Warfarin"] = "Intermediate Metabolizer"
                    raw.at[gene_index, "Warfarin Hap1"] = variant1
                else:
                    raw.at[gene_index, "Warfarin Hap1"] = "nan"
                if variant2 in variant_list:
                    if variant2 in ["*2", "*3"]:
                        raw.at[gene_index, "Phenotype Warfarin"] = "Poor Function"
                    else:
                        raw.at[gene_index, "Phenotype Warfarin"] = "Intermediate Metabolizer"
                    raw.at[gene_index, "Warfarin Hap2"] = variant2
                else:
                    raw.at[gene_index, "Warfarin Hap2"] = "nan"

            if gene == "DPYD":
                dpyd = all.loc[all["GENE"] == "DPYD"]
                target1, target2 = dpyd.loc[(dpyd["Variant 1"] == variant1)&(dpyd["Variant 2"] == variant2)], dpyd.loc[(dpyd["Variant 1"] == other_names1)&(dpyd["Variant 2"] == other_names2)]
                target3, target4 = dpyd.loc[(dpyd["Variant 1"] == variant1)&(dpyd["Variant 2"] == other_names2)], dpyd.loc[(dpyd["Variant 1"] == other_names1)&(dpyd["Variant 2"] == variant2)]
                if not target1.empty:
                    raw.at[gene_index, "Phenotype"] = target1["Phenotype"].item()
                elif not target2.empty:
                    raw.at[gene_index, "Phenotype"] = target2["Phenotype"].item()
                elif not target3.empty:
                    raw.at[gene_index, "Phenotype"] = target3["Phenotype"].item()
                elif not target4.empty:
                    raw.at[gene_index, "Phenotype"] = target4["Phenotype"].item()

        # rsid check
        for gene in exceptions_gene_rsid:

            # extract gene index from raw
            gene_index = raw.loc[raw.GENE == gene].index[0]

            # variant 1 and variant 2
            rsid1, rsid2= raw.loc[raw.GENE == gene]["rsid1"].item(), raw.loc[raw.GENE == gene]["rsid2"].item()

            # extract gene variant list
            rsid_list = exceptions_rsid.loc[exceptions_rsid.GENE == gene]["rsID"].values

            if gene == "IFNL3":
                raw.at[gene_index, "Phenotype"] = "Favorable Response"
                if rsid1 in rsid_list:
                    raw.at[gene_index, "Phenotype"], raw.at[gene_index, "hap1_main"] = "Unfavorable Response", "rs12979860(T)"
                else:
                    raw.at[gene_index, "hap1_main"] = "rs12979860 reference(C)"
                if rsid2 in rsid_list:
                    raw.at[gene_index, "Phenotype"], raw.at[gene_index, "hap2_main"] = "Unfavorable Response", "rs12979860(T)"
                else:
                    raw.at[gene_index, "hap2_main"] = "rs12979860 reference(C)"
            if gene == "VKORC1":
                raw.at[gene_index, "Phenotype"] = "Normal Metabolizer"
                if rsid1 in rsid_list:
                    raw.at[gene_index, "hap1_main"], raw.at[gene_index, "Phenotype"] = "rs9923231(T)", "Decreased Function"
                else:
                    raw.at[gene_index, "hap1_main"] = "rs9923231 reference(C)"
                if rsid2 in rsid_list:
                    raw.at[gene_index, "hap2_main"], raw.at[gene_index, "Phenotype"] = "rs9923231(T)", "Decreased Function"
                else:
                    raw.at[gene_index, "hap2_main"] = "rs9923231 reference(C)"

            if gene == "CYP4F2":
                raw.at[gene_index, "Phenotype"] = "Normal Metabolizer"
                if rsid1 in rsid_list:
                    raw.at[gene_index, "hap1_main"], raw.at[gene_index, "Phenotype"] = "rs2108622(T)","Increased Function"
                else:
                    raw.at[gene_index, "hap1_main"] = "rs2108622 reference(C)"
                if rsid2 in rsid_list:
                    raw.at[gene_index, "hap2_main"], raw.at[gene_index, "Phenotype"] = "rs2108622(T)", "Increased Function"
                else:
                    raw.at[gene_index, "hap2_main"] = "rs2108622 reference(C)"

        # fill Phenotypes with unknowns
        raw[["Phenotype"]] = raw[["Phenotype"]].fillna("Unknown")
        self.raw = raw
        return self

    def __get_pharmgkb_recommendations(self):

        self.__add_other_names()
        self.__add_rsid()
        self.__get_pharmgkb_phenotypes()

        recom = recommendations
        raw = self.raw
        all_identifier, all_medications, all_phenotypes = [], [], []
        warfarin_identifier, warfarin_medications, warfarin_phenotypes = [], [], []
        final_identifier = []

        for index, row in raw.iterrows():
            # retrieve identifying information
            gene = row["GENE"]
            variant1, variant2 = row["hap1_main"], row["hap2_main"]
            rsid1, rsid2 = row["rsid1"], row["rsid2"]
            other_names1, other_names2 = row["other_names1"], row["other_names2"]
            identifier = []

            if gene in other_names:
                if other_names1 != ".":
                    variant1 = other_names1
                if other_names2 != ".":
                    variant2 = other_names2

            genotype = variant1 + "/" + variant2
            identifier = [gene, genotype, rsid1, rsid2]
            identifier = [x for x in identifier if x != "" and not isinstance(x,float)]
            identifier = "\n".join(identifier)

            # retrieve medications
            gene_recom = recom.loc[recom.Gene == gene]
            gene_med = gene_recom["Medication(s)"].unique()
            gene_med = "\n".join(gene_med)

            # retrieve phenotype
            gene_phenotype = row["Phenotype"]

            if gene in warfarin_genes:
                if gene == "CYP2C9":
                    all_identifier.append(identifier)
                    all_phenotypes.append(gene_phenotype)
                    all_medications.append(gene_med.replace("\nWarfarin", ""))
                    gene_phenotype = row["Phenotype Warfarin"]

                warfarin_identifier.append(identifier)
                warfarin_medications.append("Warfarin")
                warfarin_phenotypes.append(gene_phenotype)
            else:
                all_identifier.append(identifier)
                all_phenotypes.append(gene_phenotype)
                all_medications.append(gene_med)

            final_identifier.append(identifier)

        # check warfarin rsID
        rsid_warfarin_table, rsid_variants, rsid_phenotype, rsid_gene, rsid_warfarin = self.rsid, ["G", "G"], "Normal Metabolizer", "RSID", "rs12777823"

        if not rsid_warfarin_table.empty:
            for index, row in rsid_warfarin_table.iterrows():
                alt = row[4]
                rsid_variants[index] = alt
                if alt == "A":
                    rsid_phenotype = "Poor Function"
                else:
                    rsid_phenotype = "Unknown"
        rsid_genotype = rsid_variants[0] + "/" + rsid_variants[1]
        rsid_identifier = [rsid_warfarin, rsid_genotype]
        rsid_identifier = "\n".join(rsid_identifier)
        rsid_warfarin_df = pd.DataFrame({"GENE": [rsid_gene], "Identifier": [rsid_identifier], "Phenotype": [rsid_phenotype]})

        report = pd.DataFrame({"Gene\nGenotype\nrsID": all_identifier + warfarin_identifier + [rsid_identifier],"Medication(s)":all_medications + warfarin_medications + ["Warfarin"],
            "Metabolizer Phenotype*":all_phenotypes + warfarin_phenotypes + [rsid_phenotype]})

        raw["Identifier"] = final_identifier

        self.raw = raw
        self._report = report
        self.rsid = rsid_warfarin_df
        return self

    def get_appendix(self):
        self.__get_pharmgkb_recommendations()

        recom = recommendations
        raw = self.raw

        #Activity Score
        activity_score_list = ["CYP2C9", "DPYD", "CYP2D6"]
        activity_score = raw.loc[raw.GENE.isin(activity_score_list)]
        activity_recom = recom.loc[recom["Medication(s)"] != "Warfarin"]
        activity_score = activity_score.merge(activity_recom, left_on=["GENE", "Phenotype", "Activity Score"], \
            right_on=["Gene","Metabolizer Phenotype*", "Activity Score"])

        # to exclude from CYP2D6
        activity_score_med = activity_score.loc[activity_score.GENE == "CYP2D6"]["Medication(s)"].values

        #G6PD
        g6pd = raw.loc[raw.GENE == "G6PD"]
        g6pd_recom = recom.loc[recom["Gene"] == "G6PD"]
        g6pd_female_recom = g6pd.copy().merge(g6pd_recom, left_on=["GENE", "Phenotype Female"], right_on=["Gene", "Metabolizer Phenotype*"])
        g6pd_female_recom["Metabolizer Phenotype*"] = g6pd_female_recom["Phenotype Female"].apply(lambda x: "Female:\n" + x)
        g6pd_male_recom = g6pd.copy().merge(g6pd_recom, left_on=["GENE", "Phenotype Male"], right_on=["Gene", "Metabolizer Phenotype*"])
        g6pd_male_recom["Metabolizer Phenotype*"] = g6pd_male_recom["Phenotype Male"].apply(lambda x: "Male:\n" + x)
        g6pd = pd.concat([g6pd_female_recom, g6pd_male_recom])

        #CYP2C9 Warfarin
        cyp2c9 = raw.loc[raw.GENE== "CYP2C9"]
        cyp2c9_hap1_recom, cyp2c9_hap2_recom = cyp2c9.copy(), cyp2c9.copy()
        cyp2c9_hap1, cyp2c9_hap2, warfarin_cyp2c9_phenotype = cyp2c9["Warfarin Hap1"].item(), cyp2c9["Warfarin Hap2"].item(), cyp2c9["Phenotype Warfarin"].item()
        recom_cyp2c9_warfarin = recom.loc[(recom.Gene == "CYP2C9") & (recom["Medication(s)"] == "Warfarin")]
        if warfarin_cyp2c9_phenotype == "Normal Metabolizer":
            cyp2c9 = cyp2c9.merge(recom_cyp2c9_warfarin, left_on=["GENE", "Phenotype Warfarin"], right_on=["Gene", "Metabolizer Phenotype*"])
        else:
            if cyp2c9_hap1 != "nan":
                cyp2c9_hap1_recom = cyp2c9_hap1_recom.merge(recom_cyp2c9_warfarin, left_on=["GENE", "hap1_main"], right_on=["Gene", "Variant"])
                cyp2c9 = cyp2c9_hap1_recom
            if cyp2c9_hap2 != "nan":
                cyp2c9_hap2_recom = cyp2c9_hap2_recom.merge(recom_cyp2c9_warfarin, left_on=["GENE", "hap2_main"], right_on=["Gene", "Variant"])
                if "Variant" in cyp2c9.columns:
                    cyp2c9 = pd.concat([cyp2c9, cyp2c9_hap2_recom])
                else:
                    cyp2c9 = cyp2c9_hap2_recom

        # all else
        exceptions_list = ["CYP2C9", "DPYD", "G6PD"]
        unknown_list, unknown = raw.loc[raw["Phenotype"].isin(["Unknown", "Female:\nUnknown\nMale:\nUnknown"])]["GENE"].tolist(), raw.loc[raw["Phenotype"].isin(["Unknown", "Female:\nUnknown\nMale:\nUnknown"])]
        unknown_recom = recom.loc[(recom["Gene"].isin(unknown_list)) & (recom["Medication(s)"] != "Warfarin")][["Gene","Medication(s)"]].drop_duplicates()
        unknown = unknown.merge(unknown_recom, left_on=["GENE"], right_on=["Gene"])
        unknown["Metabolizer Phenotype*"] = "Unknown"
        unknown["Recommendation**"] = alternative_recommendations.loc[alternative_recommendations["Metabolizer Phenotype"] == "Unknown"]["Blurb"].item()

        exceptions_list = exceptions_list + unknown_list
        recom_all = recom.loc[~(recom["Medication(s)"].isin(activity_score_med))]
        all_else = raw.loc[~(raw["GENE"].isin(exceptions_list))]
        all_else = all_else.merge(recom_all, left_on=["GENE", "Phenotype"], right_on=["Gene", "Metabolizer Phenotype*"])

        # combine all
        appendix = pd.concat([activity_score, g6pd, cyp2c9, all_else, unknown])

        # concat warfarin case to end
        non_warfarin = appendix.loc[appendix["Medication(s)"] != "Warfarin"]
        warfarin = appendix.loc[appendix["Medication(s)"] == "Warfarin"]

        # check warfarin rsID
        rsID = self.rsid.copy()
        rsid_warfarin_df = rsID.merge(recom, left_on=["GENE", "Phenotype"], right_on=["Gene", "Metabolizer Phenotype*"])

        final_appendix = non_warfarin.sort_values(by=["GENE"])
        final_appendix = pd.concat([final_appendix, warfarin, rsid_warfarin_df])
        final_appendix = final_appendix[["Identifier", "Medication(s)", "Metabolizer Phenotype*", "Implications", "Recommendation**", "Strength of Recommendation", "PharmGKB ID"]]
        final_appendix["Identifier"] = final_appendix.apply(lambda x: x["Identifier"] + "\n" + x["PharmGKB ID"] if not isinstance(x["PharmGKB ID"], float) else x["Identifier"], axis=1)
        final_appendix = final_appendix.drop(columns=["PharmGKB ID"]).fillna("N/A")
        final_appendix = final_appendix.rename(columns={"Identifier":"Gene\nGenotype\nrsID\nPharmGKB ID"})

        self._appendix = final_appendix
        return self

    def get_hla(self, hla, na):
        raw = self.raw
        report = self._report
        recom = recommendations
        if na:
            self._hla = pd.read_csv(hla, sep="\t", dtype="object").fillna("")
        else:
            self._hla = pd.read_csv(hla, sep="\t", dtype="object")

        # HLA
        hla = self._hla.copy()
        hla = hla.loc[hla.Gene.isin(["A","B"])]
        hla["GENE"] = "HLA-"+ hla.Gene
        hla["hap1_main"], hla["hap2_main"] = hla["Allele1"].apply(lambda x: x[1:7]), hla["Allele2"].apply(lambda x: x[1:7])
        hla_a, hla_b = hla.loc[hla.Gene == "A"], hla.loc[hla.Gene == "B"]
        hla_variants_a = hla_a["hap1_main"].values.tolist() + hla_a["hap2_main"].values.tolist()
        hla_variants_b = hla_b["hap1_main"].values.tolist() + hla_b["hap2_main"].values.tolist()
        hla_variants_a, hla_variants_b = set(hla_variants_a), set(hla_variants_b)

        recom_hla = recom.loc[recom.Gene.isin(["HLA-A", "HLA-B"])].copy()
        recom_hla_variants_a, recom_hla_variants_b = recom_hla["Variant"].values.tolist(), recom_hla["Variant"].values.tolist()
        recom_hla_variants_a, recom_hla_variants_b = set(recom_hla_variants_a), set(recom_hla_variants_b)
        intersection_a, intersection_b = hla_variants_a.intersection(recom_hla_variants_a), hla_variants_b.intersection(recom_hla_variants_b)

        hla_a_index, hla_b_index = hla_a.index[0].item(), hla_b.index[0].item()

        if len(intersection_a) == 0:
            hla.at[hla_a_index,"Metabolizer Phenotype*"] = "Normal or Reduced Risk"
        else:
            hla.at[hla_a_index,"Metabolizer Phenotype*"] = "Increased Risk"
        if len(intersection_b) == 0:
            hla.at[hla_b_index,"Metabolizer Phenotype*"] = "Normal or Reduced Risk"
        else:
            hla.at[hla_b_index,"Metabolizer Phenotype*"] = "Increased Risk"

        all_identifier, all_medications, all_phenotypes = [], [], []

        for index, row in hla.iterrows():
            # retrieve identifying information
            gene = row["GENE"]
            genotype = row["hap1_main"] + "/" + row["hap2_main"]
            identifier = []

            identifier = [gene, genotype]
            identifier = "\n".join(identifier)
            all_identifier.append(identifier)

            # retrieve medications
            gene_recom = recom.loc[recom.Gene == gene]
            gene_med = gene_recom["Medication(s)"].unique()
            gene_med = "\n".join(gene_med)
            all_medications.append(gene_med)

            # retrieve phenotype
            gene_phenotype = row["Metabolizer Phenotype*"]
            all_phenotypes.append(gene_phenotype)

        hla_report = pd.DataFrame({"Gene\nGenotype\nrsID": all_identifier, "Medication(s)":all_medications,
            "Metabolizer Phenotype*":all_phenotypes})
        hla["Gene\nGenotype\nrsID"] = all_identifier
        non_warfarin_report, warfarin_report = report.loc[~(report["Medication(s)"] == "Warfarin")], report.loc[report["Medication(s)"] == "Warfarin"]
        report = pd.concat([non_warfarin_report,hla_report]).sort_values(by=["Gene\nGenotype\nrsID"])
        report = pd.concat([report, warfarin_report])
        self._report = report

        # HLA
        hla_a, hla_b = hla.loc[hla.Gene == "A"].copy(), hla.loc[hla.Gene == "B"].copy()
        hla_variants_a, hla_variants_b = hla_a["hap1_main"].values.tolist() + hla_a["hap2_main"].values.tolist(), hla_b["hap1_main"].values.tolist() + hla_b["hap2_main"].values.tolist()

        recom_hla = recom.loc[recom.Gene.isin(["HLA-A", "HLA-B"])].copy()
        norm_hla, increased_hla = recom_hla.loc[recom_hla["Metabolizer Phenotype*"] == "Normal or Reduced Risk"], recom_hla.loc[recom_hla["Metabolizer Phenotype*"] == "Increased Risk"]
        norm_hla_index, increased_hla_index, final_index = norm_hla.index.tolist(), increased_hla.index.tolist(), []

        # obtain final index
        for i in range(len(norm_hla_index)):
            norm, increased = norm_hla_index[i], increased_hla_index[i]
            norm_row, increased_row = recom_hla.loc[[norm]], recom_hla.loc[[increased]]
            variant, gene = norm_row["Variant"].item(), norm_row["Gene"].item()

            if gene == "HLA-A":
                if variant in hla_variants_a:
                    final_index.append(increased)
                else:
                    final_index.append(norm)
            else:
                if variant in hla_variants_b:
                    final_index.append(increased)
                else:
                    final_index.append(norm)

        final_hla = recom_hla.loc[final_index]
        hla = hla.merge(final_hla, left_on=["GENE"], right_on=["Gene"])
        hla = hla.rename(columns={"Metabolizer Phenotype*_y":"Metabolizer Phenotype*"})
        hla["Gene\nGenotype\nrsID\nPharmGKB ID"] = hla.apply(lambda x: x["Gene\nGenotype\nrsID"] + "\n" + x["PharmGKB ID"] if not isinstance(x["PharmGKB ID"], float) else x["Gene\nGenotype\nrsID"], axis=1)

        appendix = self._appendix
        warfarin, non_warfarin = appendix.loc[appendix["Medication(s)"] == "Warfarin"], appendix.loc[~(appendix["Medication(s)"] == "Warfarin")]
        final_appendix = pd.concat([non_warfarin, hla]).sort_values(by=["Gene\nGenotype\nrsID\nPharmGKB ID"]).fillna("N/A")
        final_appendix = pd.concat([final_appendix,warfarin])
        final_appendix = final_appendix[["Gene\nGenotype\nrsID\nPharmGKB ID", "Medication(s)", "Metabolizer Phenotype*", "Implications", "Recommendation**", "Strength of Recommendation"]]
        self._appendix = final_appendix
        return self
