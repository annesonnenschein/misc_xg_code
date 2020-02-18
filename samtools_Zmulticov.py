import contextlib
import glob
import numpy as np
import pandas


def uuid_dict(input_filename):
    """sample name to uuid conversion"""
    mydict = {}
    with open(input_filename) as input_file:
        input_file.readline()
        for line in input_file:
            linelist = line.split(",")
            mydict[linelist[0]] = [linelist[0], linelist[1]]
            mydict[linelist[1]] = [linelist[0], linelist[1]]
    return mydict


def get_multicov(mydict):
    """filenames for each sample id for samtools data and 15 genes multicov"""
    myuuids = list(set([mydict[i][1] for i in mydict.keys()]))
    file_dict = {}
    multicov_file_list = glob.glob("*subpanel.bedtools_coverage_mean.tsv")
    samtools_file_list = glob.glob("*samtools_flagstat.csv")
    for each in myuuids:
        for other in multicov_file_list:
            if each in other:
                file_dict[each] = [other]
        for other in samtools_file_list:
            if each in other:
                file_dict[each].append(other)
    return file_dict


def parse_multicov(mydict, file_dict):
    """make the input from get_multicov into a dictionary"""
    sample_dict = {}
    for each_id in file_dict.keys():
        sample_id = mydict[each_id][0]
        multicov_filename = file_dict[each_id][0]
        flagstat_filename = file_dict[each_id][1]
        with open(multicov_filename) as multicov, open(flagstat_filename) as flagstat:
            flagstat_input = flagstat.readlines()
            total_reads = flagstat_input[1].strip().split(",")[-1]
            multicov_input = multicov.readlines()
            mylist = [[i.strip(), total_reads] for i in multicov_input]
            sample_dict[sample_id] = mylist
    return sample_dict


def output_files(sample_dict):
    """output total read normalized zcores"""
    header_string = "chrom,start,end,gene"
    probe_list = []
    probe_dict = {}
    for eachkey in sample_dict.keys():
        header_string = header_string + "," + eachkey
        newoutfilename = eachkey + "_multicov_flagstat.tsv"
        with open(newoutfilename, "w") as quickout:
            for eachrow in sample_dict[eachkey]:
                mean_cov = float(eachrow[0].split("\t")[-1])
                probe = ",".join(eachrow[0].split("\t")[0:4])
                total_reads = float(eachrow[1])
                outstring = probe + "," + str(mean_cov) + "," + str(total_reads) + "\n"
                quickout.write(outstring)
                normalized_cov = mean_cov / total_reads
                if probe not in probe_list:
                    probe_list.append(probe)
                if probe not in probe_dict.keys():
                    probe_dict[probe] = [normalized_cov]
                else:
                    probe_dict[probe].append(normalized_cov)
    total_out = "prevalidation_Zcoverage_by_total_reads_scale.csv"
    with open(total_out, "w") as outfile:
        outfile.write(header_string + "\n")
        for each_probe in probe_list:
            norm_mean = np.mean(probe_dict[each_probe])
            norm_sd = np.std(probe_dict[each_probe])
            z_scores = [
                str(((i - norm_mean) / norm_sd)) for i in probe_dict[each_probe]
            ]
            myoutstring = each_probe + "," + ",".join(z_scores) + "\n"
            outfile.write(myoutstring)


mydict = uuid_dict("options.list.preval")
file_dict = get_multicov(mydict)
sample_dict = parse_multicov(mydict, file_dict)
output_files(sample_dict)
