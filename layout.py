#!/usr/bin/env python
#Script shared with me when I was rotating in the Rokhsar lab, probably from Jessen Bredeson

import os
import sys
import heapq
import pysam
import getopt
import subprocess
import collections

num = len

PYTHONVERSION = sys.version_info[:2]
if PYTHONVERSION < (3,0):
    range = xrange


def which(filename):
    if os.path.exists(filename):
        return filename
    for path in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(path, filename)):
                return os.path.join(path, filename)

    return None

    
    
def span_x_with_y(x, xrc, xrl, yrc, yrl):
    post = []
    # sort hits by contig size:
    for y in sorted(xrc[x][2].keys(), key=xrc[x][2].get):
        if yrc[y][1]:
            # the contig has already been placed, so
            continue
        else:
            yrc[y][1] = 1

        ylen = yrc[y][0]
        slope = xrc[x][2][y][5][0]

        if slope < 0:
            for xx in yrc[y][2]:
                ylo = yrc[y][2][xx][0][0]
                yhi = yrc[y][2][xx][1][0]
                yrc[y][2][xx][0][0] = ylen - yhi + 1
                yrc[y][2][xx][1][0] = ylen - ylo + 1
                yrc[y][2][xx][5][0] *= -1
                
        yrl.append( [y, ylen, slope, x] )

        if ylen > xrc[x][0]:
            span_x_with_y(y, yrc, yrl, xrc, xrl)
        else:
            post.append(y)

    for y in post:
        span_x_with_y(y, yrc, yrl, xrc, xrl)



def _tile_alignments(heap, clust):
    i = 0
    idx_cluster = 0 
    num_clusters = len(clust)
    while len(heap) > 0:
        heap_ = collections.deque()
        for align in heap:
            if idx_cluster >= num_clusters:
                clust.append(align[:])
                num_clusters +=1
                continue

            if align[9] == clust[idx_cluster][9]:  # same strand
                # dQ = qry inter-hsp distance
                # dR = ref inter-hsp distance
                if align[9] < 0:
                    dQ = align[3] - clust[idx_cluster][2]
                    dR = clust[idx_cluster][6] - align[7]
                    i = 0
                else:
                    dQ = align[2] - clust[idx_cluster][3] 
                    dR = align[6] - clust[idx_cluster][7]
                    i = 1

                # sys.stderr.write("%s %s %d %d %d %d\n" % (align[0], align[4], align[9], dQ, dR, abs(dR-dQ)))

                if pow(dR**2 + dQ**2, 0.5) < 1000:
                    clust[idx_cluster][2+i] = align[2+i]
                    clust[idx_cluster][6+i] = align[6+i]
                    clust[idx_cluster][8] += (align[7] - align[6] + 1)
                else:
                    heap_.append(align)
            else:
                heap_.append(align)

        idx_cluster += 1
        heap = heap_



def tile_alignments(coords):
    def _by_query_position(query):
        beg = query[3] if query[3] < query[2] else query[2]
        end = query[2] if query[3] < query[2] else query[3]
        return (query[0], beg, end)

    heap = collections.deque()
    clusters = collections.deque()
    

    prev_qry_id = None
    for align in sorted(coords, key=_by_query_position):
        qry_id = align[0]
        ref_id = align[4]
            
        if qry_id != prev_qry_id and prev_qry_id is not None:
            _tile_alignments(heap, clusters)
            # sys.stderr.write("# %s\n" % str(clusters))
            heap = collections.deque()
        
        heap.append(align)
        prev_qry_id = qry_id

    if len(heap):
        _tile_alignments(heap, clusters)
        # sys.stderr.write("# %s\n" % str(clusters))

    return clusters



def filter_alignments(coords):
    # ensure a 1 query, 1 ref policy:
    align_max = {}
    align_all = {}
    for align in coords:
        qry_id = align[0]
        ref_id = align[4]
        ref_beg = align[6]
        ref_end = align[7]
        if qry_id not in align_all:
            align_max[qry_id] = [-1, None]
            align_all[qry_id] = {}

        if ref_id not in align_all[qry_id]:
            align_all[qry_id][ref_id] = 0
        
        align_all[qry_id][ref_id] += (ref_end - ref_beg + 1)
        if align_max[qry_id][0]  < align_all[qry_id][ref_id]:
            align_max[qry_id][0] = align_all[qry_id][ref_id]
            align_max[qry_id][1] = ref_id

    filtered = []
    for align in coords:
        qry_id = align[0]
        ref_id = align[4]
        if align_max[qry_id][1] == ref_id:
            filtered.append(align)

    return filtered



def layout(coords, ref_list, qry_list, heaviest=False):
    ref_chains = {}
    qry_chains = {}
    ref_layout = []
    qry_layout = []
    for align in coords:
        qry_id  = align[0]
        qry_len = align[1]
        qry_beg = align[2]
        qry_end = align[3]

        ref_id  = align[4]
        ref_len = align[5]
        ref_beg = align[6]
        ref_end = align[7]

        aln_len = align[8]

        if ref_id not in ref_list or qry_id not in qry_list:
            continue

        ref_ori = -1 if ref_beg >  ref_end else 1
        qry_ori = -1 if qry_beg >  qry_end else 1
        slope   = -1 if ref_ori != qry_ori else 1

        ref_min = ref_beg if ref_ori == 1 else ref_end
        ref_max = ref_end if ref_ori == 1 else ref_beg

        qry_min = qry_beg if qry_ori == 1 else qry_end
        qry_max = qry_end if qry_ori == 1 else qry_beg

        if heaviest:
            if qry_id in qry_chains:
                old_ref_id = sorted(qry_chains[qry_id][2], key=lambda r: -qry_chains[qry_id][2][r][6][0])[0]
                old_val = qry_chains[qry_id][2][old_ref_id]

                if old_val[3][0] - old_val[2][0] > ref_max - ref_min:
                    continue
                else:
                    del(ref_chains[old_ref_id][2][qry_id])
                    del(qry_chains[qry_id])
        
        # the order of the sequence length and the "placed" flag element
        # positions are swapped relative to the original mummerplot code
        # in order to allow generic sorting by sequence length
        if ref_id not in ref_chains: ref_chains[ref_id] = [ref_len, 0, {}]
        if qry_id not in qry_chains: qry_chains[qry_id] = [qry_len, 0, {}]
        if qry_id not in ref_chains[ref_id][2]: ref_chains[ref_id][2][qry_id] = [ [0] ] * 7
        if ref_id not in qry_chains[qry_id][2]: qry_chains[qry_id][2][ref_id] = [ [0] ] * 7
        
        if aln_len > ref_chains[ref_id][2][qry_id][6][0]:
            ref_chains[ref_id][2][qry_id][0] = qry_chains[qry_id][2][ref_id][2] = [ ref_min ]
            ref_chains[ref_id][2][qry_id][1] = qry_chains[qry_id][2][ref_id][3] = [ ref_max ]
            ref_chains[ref_id][2][qry_id][2] = qry_chains[qry_id][2][ref_id][0] = [ qry_min ]
            ref_chains[ref_id][2][qry_id][3] = qry_chains[qry_id][2][ref_id][1] = [ qry_max ]
            ref_chains[ref_id][2][qry_id][4] = qry_chains[qry_id][2][ref_id][4] = [ slope ]
            ref_chains[ref_id][2][qry_id][5] = qry_chains[qry_id][2][ref_id][5] = [ slope ]
            ref_chains[ref_id][2][qry_id][6] = qry_chains[qry_id][2][ref_id][6] = [ aln_len ]

    
    for ref_id in sorted(ref_chains, key=ref_chains.get, reverse=True):
        span_x_with_y(ref_id, ref_chains, ref_layout, qry_chains, qry_layout)

    for ref_id in ref_list:
        ref_list[ref_id][0] = -1

    for qry_id in qry_list:
        qry_list[qry_id][0] = -1

    roff = 0
    for r in ref_layout:
        ref_id = r[0]
        ref_list[ref_id][0] = roff
        ref_list[ref_id][2] = r[1]
        roff += ref_list[ref_id][1] - 1

    for ref_id in ref_list:
        if ref_list[ref_id][0] < 0:
            ref_list[ref_id][0] = roff
            roff += ref_list[ref_id][1] - 1
            
    qoff = 0
    for q in qry_layout:
        qry_id = q[0]
        qry_list[qry_id][0] = qoff
        qry_list[qry_id][2] = q[1]
        qoff += qry_list[qry_id][1] - 1

    for qry_id in qry_list:
        if qry_list[qry_id][0] < 0:
            qry_list[qry_id][0] = qoff
            qoff += qry_list[qry_id][1] - 1

    return(ref_chains, ref_layout, qry_chains, qry_layout)



def parse_delta_header(fn):
    fh = open(fn, 'r')

    rfile = None
    qfile = None
    mummer_type = None
    valid_mummer_types = set(('NUCMER','PROMER'))
    for line in fh:
        if line.isspace(): continue

        line = line.rstrip().split(" ")
        if len(line) == 2:
            rfile = line[0]
            qfile = line[1]
        elif line[0].upper() in valid_mummer_types:
            mummer_type = line[0].upper()
            break
        else:
            raise IOError("malformed .delta file")
        
    fh.close()
            
    if rfile is None or qfile is None or mummer_type is None:
        raise IOError("malformed .delta file")

    return(mummer_type, rfile, qfile)



def get_reference_list(references, lengths):
    ref_list = collections.OrderedDict()
    num_references = num(references)
    for i in range(num_references):
        ref_list[references[i]] = [-1, lengths[i], 1]

    return ref_list


    
def usage(message=None, status=1):
    message = '' if message is None else "\nError: %s\n\n" % message 
    sys.stderr.write("Usage: %s [-1c] [--one-to-one] [--cluster] [--out <file>] <in.delta>\n%s" % (
        os.path.basename(sys.argv[0]), message))
    sys.exit(1)

    

def main(argv):
    try:
        options, arguments = getopt.getopt(argv, 'Hhc1o:', ('help','out=','one-to-one','cluster', 'heaviest')) 
    except getopt.GetoptError as message:
        usage(message)

    strict = False
    outfile = None
    cluster = False
    heaviest = False
    for flag, value in options:
        if   flag in {'-h','--help'}: usage()
        elif flag in {'-o','--out'}: outfile = value
        elif flag in {'-1','--one-to-one'}: strict = True
        elif flag in {'-c','--cluster'}: cluster = True
        elif flag in {'-H','--heaviest'}: heaviest = True

    if len(arguments) != 1:
        usage('Unexpected number of arguments')
        
    showcoords = which('show-coords')
    if showcoords is None:
        sys.stderr.write("Error: Mummer suite not found in PATH env var\n")
        sys.exit(1)

    mum_type, ref_file, qry_file = parse_delta_header(arguments[0])

    ref_fasta = pysam.FastaFile(ref_file)
    qry_fasta = pysam.FastaFile(qry_file)

    ref_list = get_reference_list(ref_fasta.references, ref_fasta.lengths)
    qry_list = get_reference_list(qry_fasta.references, qry_fasta.lengths)
    out_file = sys.stdout if outfile is None else open(outfile, 'w')

    process = subprocess.Popen("%s -rlHT %s" % (showcoords, arguments[0]), shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    
    coords = []
    for record in stdout.decode('utf-8').splitlines():
        record = record.rstrip().split("\t")

        if len(record) != 11:
            raise IOError("Could not read show-coords pipe, invalid btab format")
        
        coords.append([
            record[10],      # qry.id
            int(record[8]),  # qry.len
            int(record[2]),  # qry.beg
            int(record[3]),  # qry.end
            record[9],       # ref.id
            int(record[7]),  # ref.len
            int(record[0]),  # ref.beg
            int(record[1]),  # ref.end
            int(record[1]) - int(record[0]) + 1,  # hsp len
            -1 if int(record[2]) > int(record[3]) else 1  # slope
        ])

    if strict:
        coords = filter_alignments(coords)
    if cluster:
        coords = tile_alignments(coords)

    ref_chains, ref_layout, qry_chains, qry_layout = layout(coords, ref_list, qry_list, heaviest)
    
    for q in qry_layout:
        out_file.write("%s\t%s\t%d\t%s\n" % (q[3], q[0], q[1], '-' if q[2] < 0 else '+'))
    
    out_file.close()

    

if __name__ == '__main__':
    main(sys.argv[1:])
