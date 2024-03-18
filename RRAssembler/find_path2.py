import networkx as nx
import time
import array
import best_path



def optimal_path(segments, contig_length, max_insertion=1500, min_aln=20, max_homology=1500, ins_cost=2.0, hom_cost=1.5,
                 inter_cost=20, intra_cost=10):
    """
    The scoring has become quite complicated.
     = Last alignment score - jump cost

    where each jump has a cost:
    jump cost = S + microhomology_length + insertion_length        # S = 10 for intra, 20 for inter

    The last score takes into account mismatched and gaps

    :param mapped:  The input dataframe with candidate alignments
    :param contig_length:   Length of the contig
    :param max_insertion: Arbitrary high number
    :param min_aln: The minimum sequence which is not part of a microhomology overlap
    :param max_homology: Arbitrary high number
    :return: A dataframe with the inferred alignments
    """

    # Start at first node then start another loop running backwards through the preceeding nodes.
    # Choose the best score out of the in edges.

    # Use a special start and end node to score the first and last alignments

    pred = array.array('i', [0] * len(segments))  # Keep track of the predecessor for the node
    node_scores = array.array('f', [0] * len(segments))

    # Deal with first score
    node_scores[0] = segments[0][3] - (segments[0][1] * ins_cost)
    pred[0] = -1

    # start from segment two because the first has been scored
    for i in xrange(1, len(segments)):
        next_chrom, next_start, next_end, next_score = segments[i]

        p = -1                                  # Bug used to be hom_cost
        best_score = next_score - (next_start * ins_cost)  # All preceeding alignments skipped!

        # Walking backwards mean the search may be terminated at some point
        for j in xrange(i-1, -1, -1):

            current_chrom, current_start, current_end, current_score = segments[j]

            # Allow alignments with minimum sequence and max overlap
            if next_start > current_end - max_homology and next_end > current_end + min_aln and \
                                    next_start - current_start > min_aln:

                if next_start > current_end and next_start - current_end > max_insertion:
                    continue

                micro_h = current_end - next_start
                if micro_h < 0:
                    ins = abs(micro_h)
                    micro_h = 0
                else:
                    ins = 0

                if next_chrom == current_chrom:
                    jump_cost = intra_cost
                else:
                    jump_cost = inter_cost

                # Calculate score, last_score = node_scores[j]
                current_score = node_scores[j] - (micro_h * hom_cost) - (ins * ins_cost) - jump_cost + next_score

                if current_score > best_score:
                    best_score = current_score
                    p = j
                    #k = (node_scores[j], (micro_h * hom_cost), (ins * ins_cost), jump_cost, next_score)

        node_scores[i] = best_score
        pred[i] = p


    # Update the score for jumping to the end of the sequence
    max_s = -1e6
    end_i = 0
    scores = []
    for i in range(len(segments)):
        sc = node_scores[i] - (ins_cost * (contig_length - segments[i][2]))
        scores.append(sc)
        # node_scores[i] = node_scores[i] - (ins_cost * (contig_length - segments[i][2]))
        if sc > max_s:
            max_s = sc
            end_i = i

    # Get path from max

    indexes = [end_i]  # [node_scores.index(max(node_scores))]

    while True:
        # Use len(indexes) - 1 to get around wraparound constraint
        next_i = pred[indexes[len(indexes) - 1]]
        if next_i == -1:
            break
        indexes.append(next_i)

    return indexes[::-1], max_s

    #
    # indexes = [end_i] #[node_scores.index(max(node_scores))]
    #
    # print "Indexes", indexes
    # while True:
    #     next_i = pred[indexes[-1]]
    #     if pred[indexes[-1]] == -1:
    #         break
    #     indexes.append(next_i)
    #
    # print max_s
    # print node_scores
    # print scores
    # scores = [node_scores[i] for i in indexes[::-1]]
    # print "Scores",  max_s / contig_length
    #
    # return [segments[i] for i in indexes[::-1]]





def old(segments, contig_length, max_insertion=1500, min_aln=20, max_homology=1500):
    G = nx.DiGraph()

    start_nodes = set([])
    end_nodes = set([])
    middle_nodes = set([])

    for i in range(len(segments)):
        u = current_chrom, current_start, current_end, current_score = segments[i]

        count = 0
        for j in range(i + 1, len(segments)):
            v = next_chrom, next_start, next_end, next_score = segments[j]

            # print current_end, next_start, next_start - current_end
            # # Skip paths with big insertions

            if next_start > current_end and next_start - current_end > max_insertion:
                continue

            # Allow alignments with minimum sequence and max overlap
            if next_start > current_end - max_homology and next_end > current_end + min_aln and \
                                    next_start - current_start > min_aln:

                # Don't connect to same sized alignments
                if (current_start, current_end) == (next_start, next_end):
                    continue

                # Add a weight parameter which describes the cost of a tranisition
                if current_chrom == next_chrom:
                    jump_cost = 10
                else:
                    jump_cost = 20

                microhomology_cost = 0
                if current_end > next_start:
                    microhomology_cost = (current_end - next_start) * 1.5

                insertion_cost = 0
                if next_start > current_end:
                    insertion_cost = (next_start - current_end) * 2.0

                total_cost = jump_cost + microhomology_cost + insertion_cost

                # Weight is the score gained by the jump (next_score) minus the cost
                # Negative and positive weights are allowed
                weight = (next_score - total_cost) * -1

                G.add_edge(u, v, weight=weight)

                if v not in middle_nodes:
                    end_nodes.add(v)

                ance_u = nx.ancestors(G, u)
                if len(ance_u) == 0:
                    start_nodes.add(u)

                else:
                    middle_nodes.add(u)
                    if u in end_nodes:
                        end_nodes.remove(u)

                count += 1

    if len(start_nodes) == 0:  # Only a single alignment or some such
        return []

    # Connect start node to each alignment. Inefficient but allows any alignment to be the start alignment.
    for n in G.nodes():
        d = (n[3] - (2 * n[1])) * -1  # Multiply by -1 because - * - is +
        G.add_edge("start", n, weight=d)

    # Connect end nodes to final end node
    for n in end_nodes:
        if n[2] == contig_length:
            distance = 0
        else:
            distance = (contig_length - n[2])
        G.add_edge(n, "end", weight=distance)

    # Shortest path algorithm on edges with negative weight is the same as the longest path problem
    pred, dist = nx.bellman_ford(G, "start", weight="weight")

    # Find the optimal path by tracing the predecessors
    p = [pred['end']]

    while True:
        parent = pred[p[-1]]
        if parent == "start":
            break
        p.append(parent)

    b = p[::-1]  # Best path
    return b


segments = (('chr1', 0, 10, 10), ('chr1', 0, 50, 50), ('chr1', 50, 109, 50), ('chr2', 55, 179, 150))

segments2 = [('chr1', 11, 61, 39), ('chr10', 44, 179, 82), ('chr14', 46, 179, 79), ('chr20', 46, 179, 96),
            ('chr6', 46, 179, 84), ('chr4', 47, 179, 90), ('chr4', 48, 179, 92), ('chr1', 49, 176, 89),
            ('chr22', 50, 176, 78), ('chr9', 50, 179, 94), ('chr6', 52, 176, 90), ('chr16', 53, 176, 87),
            ('chr19', 54, 176, 86), ('chr13', 55, 179, 81), ('chr17', 55, 179, 88), ('chr4', 55, 179, 84),
            ('chr14', 56, 179, 86), ('chr17', 56, 179, 88), ('chr19', 57, 179, 87), ('chr12', 70, 179, 82),
            ('chr10', 71, 179, 81), ('chr20', 71, 179, 80), ('chr11', 72, 179, 80), ('chr6', 72, 179, 81),
            ('chr8', 72, 176, 82)]

segments2 = sorted(segments2, key=lambda x: x[1])

indexes = (0, 18)

#print segments[indexes[0]], segments[indexes[1]]


if __name__ == "__main__":
    print

    import numpy as np
    # Change chromosome names to ints
    keys = {k: i for i, k in enumerate(set([i[0] for i in segments2]))}
    segs = np.array([[keys[k[0]], k[1], k[2], k[3]] for k in segments2])

    neg = [[0, 0, 54, 41], [1, 10, 1381, 1337], [2, 90, 1381, 1116], [3, 91, 1381, 1127], [4, 915, 1381, 399]]
    neg = [tuple(i) for i in neg]
    neg = sorted(neg, key=lambda x: x[1])
    segs = np.array(neg)
    contig_l = 1381
    print("Hi")

    t0 = time.time()
    for i in range(1):
        oldr = old(neg, contig_length=contig_l)
    t1 = time.time() -t0
    print(oldr)
    print(t1)
    # print oldr

    t0 = time.time()
    for i in range(1):
        path = optimal_path(segs, contig_l)
    print(time.time() - t0)
    print(path)


    t0 = time.time()
    for i in range(1):
        path, score = best_path.optimal_path(segs, contig_l)
    t2 = time.time() - t0

    print([segs[i] for i in path])
    print(t2)
    print(path)
    print(t1 / t2)


    #
    # import line_profiler
    # profile = line_profiler.LineProfiler(optimal_path)
    # profile.runcall(optimal_path, segments2, 179)
    # profile.print_stats()
