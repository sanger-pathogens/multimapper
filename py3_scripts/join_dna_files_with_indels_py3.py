#!/usr/bin/env python

from optparse import OptionParser
from modules.Si_SeqIO import *


#################################
# Simple Error Printing Funtion #
#################################

def DoError(ErrorString):
    print("!!!Error:", ErrorString, "!!!")
    sys.exit()


##########################################
# Function to Get command line arguments #
##########################################


def main():
    usage = "usage: %prog [options] <list of mfa files>"
    parser = OptionParser(usage=usage)

    parser.add_option("-r", "--reference", action="store", dest="ref",
                      help="Reference fasta (or multifasta). Must be the one used for mapping", default="")

    parser.add_option("-o", "--output", action="store", dest="output", help="output file prefix", default="")
    parser.add_option("-c", "--curate", action="store_true", dest="curate",
                      help="Manually curate added insertions using seaview [default=%default]", default=False)
    parser.add_option("-t", "--textfile", action="store", dest="textfile", help="textfile containing input mfa files",
                      default="")

    return parser.parse_args()


################
# Main program #
################		

if __name__ == "__main__":

    # Get command line arguments

    (options, args) = main()

    if options.ref == "":
        DoError("Reference file (-r) required")

    if options.output == "":
        DoError("Output file (-o) required")

    if options.textfile == "":
        DoError("No textfile of input mfa files selected")
    if not os.path.isfile(options.textfile):
        DoError("Cannot find file " + options.textfile)

    # Read the reference file

    if not os.path.isfile(options.ref):
        DoError("Cannot find file " + options.ref)

    # read the reference fasta file

    try:
        fasta = SeqIO.parse(open(options.ref), "fasta")
    except Exception:
        DoError("Cannot open file " + options.ref)

    if options.textfile != "" and not os.path.isfile(options.textfile):
        DoError("Cannot find file " + options.textfile)

    for line in open(options.textfile, "rU"):
        line = line.strip()
        args.append(line)

    reforder = []
    refseq = {}
    for sequence in fasta:
        if sequence.id not in refseq:
            refseq[sequence.id] = {}
            reforder.append(sequence.id)
        refseq[sequence.id] = sequence
        refseq[sequence.id].id = '.'.join(options.ref.split("/")[-1].split(".")[:-1])

    # Get all insertion locations
    Insertions = {}
    Insertion_locations = {}
    Deletions = {}
    Deletion_locations = {}

    for arg in args:
        if arg.split(".")[-1] == "mfa":
            indelfile = arg.replace(".mfa", "_indels.txt")
        elif arg.split(".")[-1] == "dna":
            indelfile = arg.replace(".dna", "_indels.txt")
        else:
            DoError("Input files must end in .dna or .mfa")
        if not os.path.isfile(indelfile):
            print("Cannot find ", indelfile)
        else:
            try:
                lines = open(indelfile, "rU").readlines()
            except Exception:
                print("Cannot open ", indelfile)

            for line in lines:
                contig = line.split()[0]
                location = int(line.split()[1]) - 1
                # note that insertions must be +1, as they are AFTER the base
                indeltype = line.strip().split()[2]
                change = line.strip().split()[3]
                if indeltype == "+":
                    if contig not in Insertion_locations:
                        Insertion_locations[contig] = []
                        Insertions[contig] = {}

                    if location not in Insertions[contig]:
                        Insertion_locations[contig].append(location)
                        Insertions[contig][location] = {}
                    Insertions[contig][location][arg] = change
                elif indeltype == "-":
                    if contig not in Deletion_locations:
                        Deletion_locations[contig] = []
                        Deletions[contig] = {}

                    if location not in Deletions[contig]:
                        Deletion_locations[contig].append(location)
                        Deletions[contig][location] = {}
                    Deletions[contig][location][arg] = change

    # Read the mfa files for each isolate
    sequences = {}

    for arg in args:

        try:
            fasta = SeqIO.parse(open(arg), "fasta")
        except Exception:
            DoError("Cannot open file " + arg)

        for sequence in fasta:

            seqid = ""

            if sequence.id in reforder:
                seqid = sequence.id
            else:
                x = 0
                while x < len(sequence.id.split("_")):

                    if "_".join(sequence.id.split("_")[x:]) in reforder:
                        seqid = "_".join(sequence.id.split("_")[x:])
                        break
                    # print reforder, sequence.id,  seqid
                    x += 1
            if seqid == "":
                DoError(seqid + " and " + sequence.id + " not in reference genome")
            if seqid not in sequences:
                sequences[seqid] = {}
            sequences[seqid][arg] = sequence
            sequences[seqid][arg].id = '.'.join(arg.split("/")[-1].split(".")[:-1])

    chars = string.ascii_letters + string.digits
    tmpname = 'tmp' + "".join(choice(chars) for x in range(randint(8, 10)))

    # Deal with the deletions

    align_size = 20

    for contig in Deletion_locations:
        Deletion_locations[contig].sort()
        Deletion_locations[contig].reverse()

    for contig in sequences:
        if not contig in Deletion_locations:
            continue

        for location in Deletion_locations[contig]:
            # print location
            maxdellen = 0
            # lens=[]
            for sequence in sequences[contig]:
                if sequence in Deletions[contig][location]:
                    if len(Deletions[contig][location][sequence]) > maxdellen:
                        maxdellen = len(Deletions[contig][location][sequence])
                    # lens.append(Deletions[contig][location][sequence])

            alignlen = maxdellen + align_size

            if location >= align_size and location <= len(refseq[contig]) - (align_size + 1):
                start = location - align_size
                end = location + alignlen
            elif location < align_size:
                start = 0
                end = location + alignlen
            else:
                end = len(refseq[contig])
                start = location - align_size

            tempseqs = [refseq[contig][start:end]]

            for sequence in sequences[contig]:

                tempseqs.append(sequences[contig][sequence][start:end])
                tempseqs[-1].id = sequence
                if sequence in Deletions[contig][location]:

                    deletionlen = len(Deletions[contig][location][sequence])
                    if start == 0:
                        distfromstart = align_size - location

                    else:
                        distfromstart = align_size

                    startbit = tempseqs[-1].seq[:distfromstart + 1].upper()

                    if distfromstart + deletionlen < len(tempseqs[-1].seq):
                        endbit = tempseqs[-1].seq[distfromstart + 1 + deletionlen:].upper()
                    else:
                        endbit = ""

                    tempseqs[-1].seq = startbit + "-" * deletionlen + endbit

            SeqIO.write(tempseqs, open(tmpname + ".aln", "w"), "fasta")
            if options.curate:
                os.system("seaview " + tmpname + ".aln")

            muscleout = open(tmpname + ".aln", "rU").read().split(">")[1:]
            for line in muscleout:
                words = line.split("\n")
                name = words[0].split()[0]
                if name not in sequences[contig]:
                    refseq[contig].seq = refseq[contig].seq[:start] + ''.join(words[1:]) + refseq[contig].seq[end:]
                    continue
                else:
                    sequences[contig][name].seq = sequences[contig][name].seq[:start] + ''.join(words[1:]) + \
                                                  sequences[contig][name].seq[end:]

            refseqnewlen = len(refseq[contig].seq)

            for name in sequences[contig]:
                if len(sequences[contig][name].seq) != refseqnewlen:
                    print("D", start, end, location, refseqnewlen, len(sequences[contig][name].seq))
                    sys.exit()
                    os.system("seaview " + tmpname + ".aln")

    # sort the inserts

    for contig in Insertion_locations:
        Insertion_locations[contig].sort()
        Insertion_locations[contig].reverse()

    # Deal with the inserts (most complicated bit)

    for contig in sequences:
        if not contig in Insertion_locations:
            continue

        for location in Insertion_locations[contig]:

            if location >= align_size and location <= len(refseq[contig]) - (align_size + 1):
                start = location - align_size
                end = location + align_size
            elif location < align_size:
                start = 0
                end = align_size * 2
            else:
                end = len(refseq[contig]) - 1
                start = len(refseq[contig]) - (align_size * 2) + 1
            tempseqs = [refseq[contig][start:end]]

            for sequence in list(sequences[contig].keys()):
                tempseqs.append(sequences[contig][sequence][start:end])
                tempseqs[-1].id = sequence
                if sequence in Insertions[contig][location]:
                    tempseqs[-1].seq = tempseqs[-1].seq[:align_size + 1] + Insertions[contig][location][sequence] + \
                                       tempseqs[-1].seq[align_size + 1:]

            SeqIO.write(tempseqs, open(tmpname + ".fasta", "w"), "fasta")
            os.system("muscle -in " + tmpname + ".fasta -out " + tmpname + ".aln  >  /dev/null 2>&1 ")

            if options.curate:
                os.system("seaview " + tmpname + ".aln")
            muscleout = open(tmpname + ".aln", "rU").read().split(">")[1:]
            for line in muscleout:
                words = line.split("\n")
                name = words[0].split()[0]
                if name not in sequences[contig]:
                    refseq[contig].seq = refseq[contig].seq[:start] + ''.join(words[1:]) + refseq[contig].seq[end:]
                else:
                    seqline = ''.join(words[1:])
                    while ("N-" in seqline) or ("-N" in seqline):
                        seqline = seqline.replace("N-", "NN").replace("-N", "NN")
                    sequences[contig][name].seq = sequences[contig][name].seq[:start] + seqline + sequences[contig][
                                                                                                      name].seq[end:]

            refseqnewlen = len(refseq[contig].seq)

            for name in sequences[contig]:
                if len(sequences[contig][name].seq) != refseqnewlen:
                    print("I", start, end, location)
                    print(tempseqs)
                    sys.exit()

    final_sequences = []
    final_sequences.append(refseq[reforder[0]])
    for contig in reforder[1:]:
        final_sequences[-1].seq = final_sequences[-1].seq + refseq[contig].seq

    for sequence in list(sequences[contig].keys()):
        final_sequences.append(sequences[reforder[0]][sequence])
        for contig in reforder[1:]:
            final_sequences[-1].seq = final_sequences[-1].seq + sequences[contig][sequence].seq

    SeqIO.write(final_sequences, open(options.output, "w"), "fasta")

    os.system("rm -f " + tmpname + "*")
