# Author Gaurav
# Univeristat Potsdam
# Date 2024-5-6
# a streamlit application for the pacbiohifi from the sequencing to the read

import streamlit as st
import pandas as pd
import streamlit.components.v1 as components

st.set_page_config(
                 page_title="PacbioHifi Read Analyzer",
                 layout="centered",
                 initial_sidebar_state="expanded")
st.header("PacBioHifi Analyzer")
st.subheader("Developed by Gaurav Sablok")

help = st.button("Display the help toggle button")
if help:
    st.write("The following options are present in the Streamlit PacBioHifi application")
    st.write("1. FASTQ reader: reads pacbiohifi reads and plots them")
    st.write("2. FASTQ to FASTA converter: read the pacbiohifi reads and converts them to the fasta")
    st.write("3. FASTQ filter: filter the pacbiohifi reads with the specific clip sequences")
    st.write("4. FASTQ/FASTA length plotter")
    st.write("5. ReadChecker: checks the reads for the dna strings and plot the \
              before and after those reads, supports the output fasta file writing ")
    st.write("6. ReadExtractor: extracts the reads with the specific dna patterns and plots the before \
        and after them ")

# multi option display menu from the read analysis to the read string checker and the read extractor.

options = st.selectbox("Please select the options based on the FASTA/FASTQ",["FASTA/FASTQ Read Analysis", "FASTA/FASTQ read string checker", "FASTA/FASTQ ReadExtractor"])
if options == "FASTA/FASTQ Read Analysis":
    filetype = st.selectbox("Please select the type of the files: fastq or the fasta", ["fastq", "fasta"])
    if filetype == "fasta":
        filepath = st.text_input("enter the file path")
        read_transcripts = [i.strip() for i in open(filepath, "r").readlines()]
        fasta_dict = {}
        for i in read_transcripts:
            if i.startswith(">"):
                path = i.strip()
                if i not in fasta_dict:
                    fasta_dict[i] = ""
                    continue
            fasta_dict[path] += i.strip()
        fasta_seq = list(fasta_dict.values())
        fasta_names = [i.replace(">", "")for i in (list(fasta_dict.keys()))]
        lenfasta = []
        for i in range(len(fasta_seq)):
            lenfasta.append(len(fasta_seq[i]))
        lendata = pd.DataFrame(lenfasta, columns=["PacbioHifi Length"])
        names = st.checkbox("Press if the fasta names are needed")
        sequences = st.checkbox("Press if the fasta sequences are needed")
        length = st.checkbox("Press if the length plot are needed")
        if names:
            st.write(f"the names are: {fasta_names}")
        if sequences:
            st.write(f"the sequences are:{fasta_seq}")
        if length:
            st.bar_chart(lendata)
    if filetype == "fastq":
        filepath = st.text_input("enter the file path")
        option = st.selectbox("Please select the option where to display or write the fasta", ["display", "write"])
        if filepath and option == "display":
            readfastq = [i.strip() for i in open(filepath, "r").readlines()]
            fastq_dict = {}
            for i in range(len(readfastq)):
                if readfastq[i].startswith("@"):
                    fastq_dict[readfastq[i]] = readfastq[i + 1]
            fastq_names = list(fastq_dict.keys())
            fastq_sequences = list(fastq_dict.values())
            fastq_length = list(map(lambda n: len(n), fastq_sequences))
            length = pd.DataFrame(fastq_length, columns=["PacbioHifi Length"])
            names = st.checkbox("Press if the fasta names are needed")
            sequences = st.checkbox("Press if the fasta sequences are needed")
            lengthplot = st.checkbox("Press if the length plot are needed")
            if names:
                st.write(f"the names are: {fastq_names}")
            if sequences:
                st.write(f"the sequences are:{fastq_sequences}")
            if lengthplot:
                st.bar_chart(length)
        if filepath and option == "write":
            fileoutput = st.text_input("enter the path for the output file")
            readfastq = [i.strip() for i in open(filepath, "r").readlines()]
            fastq_dict = {}
            for i in range(len(readfastq)):
                if readfastq[i].startswith("@"):
                    fastq_dict[readfastq[i]] = readfastq[i + 1]
            fastq_names = list(fastq_dict.keys())
            fastq_sequences = list(fastq_dict.values())
            fastq_length = list(map(lambda n: len(n), fastq_sequences))
            length = pd.DataFrame(fastq_length, columns=["PacbioHifi Length"])
            names = st.checkbox("Press if the fasta names are needed")
            sequences = st.checkbox("Press if the fasta sequences are needed")
            lengthplot = st.checkbox("Press if the length plot are needed")
            if names:
                st.write(f"the names are: {fastq_names}")
            if sequences:
                st.write(f"the sequences are:{fastq_sequences}")
            if lengthplot:
                st.bar_chart(length)
            with open(fileoutput, "w") as writefile:
                for i in range(len(fastq_names)):
                    writefile.write(f">{fastq_names[i]}\n{fastq_sequences[i]}\n")

if options == "FASTA/FASTQ read string checker":
    store = st.text_input("Please enter the dna string that you want to check in the reads")
    if store and st.button("check sequence pattern"):
        st.markdown("This option is only available for the fasta files")
        readpatternplot = st.checkbox("Do you want to plot the reads with the pattern")
        if not readplot:
            filepath = st.text_input("Please enter the path for the fasta files")
            read_transcripts = [i.strip() for i in open(filepath, "r").readlines()]
            fasta_dict = {}
            for i in read_transcripts:
                if i.startswith(">"):
                    path = i.strip()
                    if i not in fasta_dict:
                        fasta_dict[i] = ""
                    continue
                fasta_dict[path] += i.strip()
            fasta_seq = list(fasta_dict.values())
            fasta_names = [i.replace(">", "")for i in (list(fasta_dict.keys()))]
            selectedones = {}
            for i in range(len(fasta_seq)):
                if store in fasta_seq[i]:
                    selectedones[fasta_names[i]] = fasta_seq[i]
            with open(filepath, "w") as writefile:
                for k,v in selectedones.items():
                    writefile.write(f">{k}\n{v}")
            st.write("The file has been written")
        if readpatternplot:
            filepath = st.text_input("Please enter the path for the fasta files")
            read_transcripts = [i.strip() for i in open(filepath, "r").readlines()]
            fasta_dict = {}
            for i in read_transcripts:
                if i.startswith(">"):
                    path = i.strip()
                    if i not in fasta_dict:
                        fasta_dict[i] = ""
                    continue
                fasta_dict[path] += i.strip()
            fasta_seq = list(fasta_dict.values())
            fasta_names = [i.replace(">", "")for i in (list(fasta_dict.keys()))]
            selectedones = {}
            for i in range(len(fasta_seq)):
                if store in fasta_seq[i]:
                    selectedones[fasta_names[i]] = fasta_seq[i]
            with open(filepath, "w") as writefile:
                for k,v in selectedones.items():
                    writefile.write(f">{k}\n{v}")
            st.write("The file has been written")
            names = list(selectedones.keys())
            selectedlength = list(map(lambda n: len(n),list(selectedones.values())))
            st.bar_chart(selectedlength)

sequencestring = st.text_input("Please enter the string pattern that you want to check in the reads")
if sequencestring and st.button("check sequence pattern"):
    st.markdown("This option is only available for the fasta files")
    readpatternplot = st.checkbox("Do you want to plot the reads with the pattern")
    filepath = st.text_input("Please enter the path for the fasta files")
    read_transcripts = [i.strip() for i in open(filepath, "r").readlines()]
    fasta_dict = {}
    for i in read_transcripts:
        if i.startswith(">"):
            path = i.strip()
            if i not in fasta_dict:
                fasta_dict[i] = ""
            continue
        fasta_dict[path] += i.strip()
    fasta_seq = list(fasta_dict.values())
    fasta_names = [i.replace(">", "")for i in (list(fasta_dict.keys()))]
    selectedones = {}
    for i in range(len(fasta_seq)):
        if sequencestring in fasta_seq[i]:
            selectedones[fasta_names[i]] = fasta_seq[i]
    with open(filepath, "w") as writefile:
        for k,v in selectedones.items():
            writefile.write(f">{k}\n{v}")
    st.write("The file has been written")
    if readpatternplot:
        filepath = st.text_input("Please enter the path for the fasta files")
        read_transcripts = [i.strip() for i in open(filepath, "r").readlines()]
        fasta_dict = {}
        for i in read_transcripts:
            if i.startswith(">"):
                path = i.strip()
                if i not in fasta_dict:
                    fasta_dict[i] = ""
                continue
            fasta_dict[path] += i.strip()
        fasta_seq = list(fasta_dict.values())
        fasta_names = [i.replace(">", "")for i in (list(fasta_dict.keys()))]
    selectedones = {}
    for i in range(len(fasta_seq)):
        if storingstring in fasta_seq[i]:
            selectedones[fasta_names[i]] = fasta_seq[i]
    with open(filepath, "w") as writefile:
        for k,v in selectedones.items():
            writefile.write(f">{k}\n{v}")
    st.write("The file has been written")
    names = list(selectedones.keys())
    selectedlength = list(map(lambda n: len(n),list(selectedones.values())))
    st.bar_chart(selectedlength)

 # filtering the reads according to the read length and plotting them before and after

readlength = st.text_input("Please enter the read threshold for the PacbioHifi reads:")
if readlength and st.button("PacbIoHifi Read Filter"):
    typefile = st.selectbox("Please select the type of the files: fastq or the fasta", ["fastq", "fasta"])
    option = st.selectbox("Please select the option:",["read", "write"])
    if typefile == "fasta" and option == "read":
        read_transcripts = [i.strip() for i in open(filepath, "r").readlines()]
        fasta_dict = {}
        for i in read_transcripts:
            if i.startswith(">"):
                path = i.strip()
                if i not in fasta_dict:
                    fasta_dict[i] = ""
                continue
            fasta_dict[path] += i.strip()
        fasta_seq = list(fasta_dict.values())
        fasta_names = [i.replace(">", "")for i in (list(fasta_dict.keys()))]
        lenfasta = []
        for i in range(len(fasta_seq)):
            lenfasta.append(len(fasta_seq[i]))
        fastacombined = {}
        for i in range(len(fasta_names)):
            fastacombined[fasta_names[i]] = [lenfasta[i], fasta_seq[i]]
        fastacombinedfiltered = [(k,v) for k,v in fastacombined.items() if v[0] > int(readlength)]
        fastacombinedfilteredlength = [v for k,v in fastacombined.items()]
        names = st.checkbox("Press if the fasta names are needed")
        sequences = st.checkbox("Press if the fasta sequences are needed")
        length = st.checkbox("Press if the length plot are needed")
        comparative = st.checkbox("Press if the comparative plots before and after filtering")
        if names:
            for k,v in fastacombinedfiltered.items():
                st.write(f"the names are: {k}")
        if sequences:
            for k,v in fastacombinedfiltered.items():
                st.write(f"the sequences are:{v}")
        if length:
            st.bar_chart(fastacombinedfilteredlength)
        if comparative:
            st.bar_chart(lenfasta)
            st.bar_chart(fastacombinedfilteredlength)
    if typefile == "fasta" and option == "write":
        filepath = st.text_input("enter the file path")
        read_transcripts = [i.strip() for i in open(filepath, "r").readlines()]
        fasta_dict = {}
        for i in read_transcripts:
            if i.startswith(">"):
                path = i.strip()
                if i not in fasta_dict:
                    fasta_dict[i] = ""
                continue
            fasta_dict[path] += i.strip()
        fasta_seq = list(fasta_dict.values())
        fasta_names = [i.replace(">", "")for i in (list(fasta_dict.keys()))]
        lenfasta = []
        for i in range(len(fasta_seq)):
            lenfasta.append(len(fasta_seq[i]))
        fastacombined = {}
        for i in range(len(fasta_names)):
            fastacombined[fasta_names[i]] = [lenfasta[i], fasta_seq[i]]
        fastacombinedfiltered = [(k,v) for k,v in fastacombined.items() if v[0] > int(readlength)]
        fastacombinedfilteredlength = [v for k,v in fastacombined.items()]
        names = st.checkbox("Press if the fasta names to be written")
        sequences = st.checkbox("Press if the fasta sequences to be written")
        length = st.checkbox("Press if the length plot are needed")
        comparative = st.checkbox("Press if the comparative plots before and after filtering")
        if names:
             with open(filewrite, "w") as writefile:
                   for k,v in fastacombinedfiltered.items():
                       writefile.write(f">{k}\n")
        if sequences:
           with open(filewrite, "w") as writefile:
                   for k,v in fastacombinedfiltered.items():
                       writefile.write(f">{v}\n")
        if length:
            st.bar_chart(fastacombinedfilteredlength)
        if comparative:
            st.bar_chart(lenfasta)
            st.bar_chart(fastacombinedfilteredlength)
        if fasta:
            with open(filewrite, "w") as writefile:
                   for k,v in fastacombinedfiltered.items():
                       writefile.write(f">{k}\n{v}\n")

    if typefile == "fastq" and option == "read":
        read_transcripts = [i.strip() for i in open(filepath, "r").readlines()]
        fasta_dict = {}
        readfastq = [i.strip() for i in open(filepath, "r").readlines()]
        fastq_dict = {}
        for i in range(len(readfastq)):
            if readfastq[i].startswith("@"):
                fastq_dict[readfastq[i]] = readfastq[i + 1]
        fastq_names = list(fastq_dict.keys())
        fastq_sequences = list(fastq_dict.values())
        fastq_length = list(map(lambda n: len(n), fastq_sequences))
        fastqcombined = {}
        for i in range(len(fastq_names)):
            fastqcombined[fastq_names[i]] = [fastq_length[i], fastq_sequences[i]]
        fastqcombinedfiltered = [(k,v) for k,v in fastacombined.items() if v[0] > int(readlength)]
        fastqcombinedfilteredlength = [v for k,v in fastqcombined.items()]
        names = st.checkbox("Press if the fasta names are needed")
        sequences = st.checkbox("Press if the fasta sequences are needed")
        length = st.checkbox("Press if the length plot are needed")
        comparative = st.checkbox("Press if the comparative plots before and after filtering")
        if names:
            for k,v in fastqcombinedfiltered.items():
                st.write(f"the names are: {k}")
        if sequences:
            for k,v in fastqcombinedfiltered.items():
                st.write(f"the sequences are:{v}")
        if length:
            st.bar_chart(fastqcombinedfilteredlength)
        if comparative:
            st.bar_chart(lenfasta)
            st.bar_chart(fastqcombinedfilteredlength)
    if typefile == "fastq" and option == "write":
        filepath = st.text_input("enter the file path")
        read_transcripts = [i.strip() for i in open(filepath, "r").readlines()]
        fasta_dict = {}
        readfastq = [i.strip() for i in open(filepath, "r").readlines()]
        fastq_dict = {}
        for i in range(len(readfastq)):
            if readfastq[i].startswith("@"):
                fastq_dict[readfastq[i]] = readfastq[i + 1]
        fastq_names = list(fastq_dict.keys())
        fastq_sequences = list(fastq_dict.values())
        fastq_length = list(map(lambda n: len(n), fastq_sequences))
        fastqcombined = {}
        for i in range(len(fastq_names)):
            fastqcombined[fastq_names[i]] = [fastq_length[i], fastq_sequences[i]]
        fastqcombinedfiltered = [(k,v) for k,v in fastacombined.items() if v[0] > int(readlength)]
        fastqcombinedfilteredlength = [v for k,v in fastqcombined.items()]
        names = st.checkbox("Press if the fasta names to be written")
        sequences = st.checkbox("Press if the fasta sequences to be written")
        length = st.checkbox("Press if the length plot are needed")
        comparative = st.checkbox("Press if the comparative plots before and after filtering")
        if names:
             with open(filewrite, "w") as writefile:
                   for k,v in fastqcombinedfiltered.items():
                       writefile.write(f">{k}\n")
        if sequences:
           with open(filewrite, "w") as writefile:
                   for k,v in fastqcombinedfiltered.items():
                       writefile.write(f">{v}\n")
        if length:
            st.bar_chart(fastqcombinedfilteredlength)
        if comparative:
            st.bar_chart(lenfasta)
            st.bar_chart(fastqcombinedfilteredlength)
        if fasta:
            with open(filewrite, "w") as writefile:
                   for k,v in fastqcombinedfiltered.items():
                       writefile.write(f">{k}\n{v}\n")

# checking patterns in the PacBioHifi Reads

patternstring = st.text_input("Please enter the string that you want to extract")
if patternstring and st.button("find the string occurences"):
    st.write("This option is available for the fasta files")
    st.write("Please convert your fastq file into fasta")
    confirmadd = st.button("Please confirm with the file notion")
    if confirmadd:
        filepath = st.text_input("Please enter the path for the fasta/fastq files")
        fileout = st.text_input("Please enter the path for the output fasta files")
        if filepath == "fasta" and fileout:
            read_transcripts = [i.strip() for i in open(filepath, "r").readlines()]
            fasta_dict = {}
            for i in read_transcripts:
                if i.startswith(">"):
                    path = i.strip()
                    if i not in fasta_dict:
                        fasta_dict[i] = ""
                    continue
                fasta_dict[path] += i.strip()
            fasta_seq = list(fasta_dict.values())
            fasta_names = [i.replace(">", "")for i in (list(fasta_dict.keys()))]
            lengthpatternstring = len(patternstring)
            fastastringfiltered = {}
            for i in range(len(fasta_names)):
                if addstring in fasta_seq[i]:
                    fastastringfiltered[fasta_names[i]] = [int(fasta_seq[i].find(patternstring)), int(fasta_seq[i].find(addstring)+lenaddstring), fasta_seq[i]]
            fastastringfilterednames = list(fastastringfiltered.keys())
            fastastringfilteredstart = [i[0] for i in fastastringfiltered.values()]
            fastastringfilteredend = [i[1] for i in fastastringfiltered.values()]
            fastastringfilteredseq = [i[2] for i in fastastringfiltered.values()]
            sliceout = {}
            for i in range(len(fastastringfilterednames)):
                sliceout[fastastringfilterednames[i]] = fastastringfilteredseq[i][fastastringfilteredstart[i]:fastastringfilteredend[i]]
                with open(fileout, "w") as writefasta:
                    for k,v in slicedout.items:
                        writefile.write(f"{k}\n{v}")
