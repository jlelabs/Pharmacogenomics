import sys
import pandas as pd
import pharma_class as pharma

report_phenotypes = ["Ultrarapid Metabolizer", "Rapid Metabolizer", "Intermediate Metabolizer",\
    "Poor Metabolizer", "Decreased Function", "Poor Function", "Increased Risk", "Unfavorable Response",\
    "Malignant Hyperthermia Susceptible", "Deficient", "Variable", "Ivacaftor responsive in CF patients",\
    "Increased Function", "Unknown"]


# testing
sample_name = sys.argv[1]
path_to_output = sys.argv[2] + "/"
path_to_sup = sys.argv[3] + "/"

path_to_stargazer = path_to_sup + sample_name + "_final_stargazer_output.txt"
rsID_warfarin = path_to_sup + sample_name + "_rsID_Warfarin.txt"

subject = pharma.SubjectPharmaBuilder(path_to_stargazer,rsID_warfarin)
subject = subject.get_appendix()

report_path = path_to_sup + "/" + sample_name + "_report.xlsx"
appendix_path = path_to_output + "/" + sample_name + "_appendix.xlsx"
subject_raw = path_to_sup + "/" + sample_name + "_subject_raw.xlsx"

subject._report.to_excel(report_path, index=False)
subject._appendix.to_excel(appendix_path, index=False)
subject.raw.to_excel(subject_raw, index=False)

non_normal_report = subject._report.copy()
non_normal_report = non_normal_report.loc[non_normal_report["Metabolizer Phenotype*"].isin(report_phenotypes)]
non_normal_report.to_excel(path_to_output + "/" + sample_name + "_report.xlsx", index=False)
