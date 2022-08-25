
#MIT License
#
#Copyright (c) 2022 LIT-University of Strasbourg
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


# Code used in
# Target-focused library design by pocket-applied computer vision and fragment deep generative linking
# Merveille Eguida, Christel Valencia, Marcel Hibert, Pascal Villa and Didier Rognan
# contacts: keguida[at]unistra.fr, rognan[at]unistra.fr



import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--ifile', type=str, required=True,
            help='input file')
parser.add_argument('-o', '--ofile', type=str, required=True,
                help='output file')
args = parser.parse_args()


with open(args.ofile, 'w') as of:
    with open(args.ifile, 'r') as f:
        i = 1
        for line in f:
            cols = line.split()
            cols.append(i)
            of.write(('\t'.join(['{}' for n in range(len(cols))])+'\n').format(*cols))
            i += 1

