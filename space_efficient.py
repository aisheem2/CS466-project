from memory_profiler import profile

from Hirschberg import hirschberg


def find_end(v, w, delta, mode="fitting"):
    assert mode in ["fitting", "local"], "The alignment must be fitting alignment or local alignment"

    M = [[0 for j in range(2)] for i in range(len(v)+1)]
    ### BEGIN SOLUTION
    def needleman(v, w, delta):
        max_score = float("-inf")
        max_score_i, max_score_j = None, None
        for j in range(len(w)+1):
            for i in range(len(v)+1):
                if i == 0:
                    M[i][j % 2] = 0
                elif j == 0:
                    if mode == "fitting":
                        M[i][j % 2] = M[i - 1][j % 2] + delta[v[i - 1]]['-']
                    elif mode == "local":
                        M[i][j % 2] = 0
                    else:
                        raise NotImplementedError()
                else:
                    best_sub = max([
                        M[i][(j - 1) % 2] + delta['-'][w[j - 1]],
                        M[i - 1][j % 2] + delta[v[i - 1]]['-'],
                        M[i - 1][(j - 1) % 2] + delta[v[i - 1]][w[j - 1]]], key = lambda x: x)
                    M[i][j % 2] = best_sub
                if mode == "fitting":
                    if (i == len(v)) and (M[i][j % 2] > max_score):
                        max_score = M[i][j % 2]
                        max_score_i, max_score_j = i, j
                elif mode == "local":
                    if M[i][j % 2] > max_score:
                        max_score = M[i][j % 2]
                        max_score_i, max_score_j = i, j
                else:
                    raise NotImplementedError()
        return max_score, max_score_i, max_score_j
    
    score, v_end, w_end = needleman(v, w, delta)
    return score, v_end - 1, w_end - 1


def find_start(v, v_end, w, w_end, delta, mode="fitting"):
    assert mode in ["fitting", "local"], "The alignment must be fitting alignment or local alignment"

    v_ = v[0 : v_end + 1]
    w_ = w[0 : w_end + 1]

    v_ = v_[::-1]
    w_ = w_[::-1]

    _, v_end, w_end = find_end(v_, w_, delta, mode)
    
    v_start = len(v_) - 1 - v_end
    w_start = len(w_) - 1 - w_end

    return v_start, w_start


@profile
def space_efficient_fitting_align(v, w, delta):
    score, v_end, w_end = find_end(v, w, delta, "fitting")
    v_start, w_start = find_start(v, v_end, w, w_end, delta, "fitting")
    aligned_v, aligned_w = hirschberg(v[v_start : v_end + 1], w[w_start : w_end + 1])

    return score, aligned_v, aligned_w


@profile
def space_efficient_local_align(v, w, delta):
    score, v_end, w_end = find_end(v, w, delta, "local")
    v_start, w_start = find_start(v, v_end, w, w_end, delta, "local")
    aligned_v, aligned_w = hirschberg(v[v_start : v_end + 1], w[w_start : w_end + 1])

    return score, aligned_v, aligned_w


if __name__ == "__main__":
    X = "TGC"
    Y = "ACTGCAA"

    keys = ['A', 'C', 'T', 'G', '-']
    delta = {}
    for i in range(len(keys)):
        delta[keys[i]] = {k : v for (k,v) in zip(keys, [1 if keys[i] == keys[j]  else -1 for j in range(len(keys))])}

    mode = "fitting"
    score, v_end, w_end = find_end(X, Y, delta, mode)
    v_start, w_start = find_start(X, v_end, Y, w_end, delta, mode)
    print("score: {}".format(score))
    print("v_start: {} v_end: {}".format(v_start, v_end))
    print("w_start: {} w_end: {}".format(w_start, w_end))

    
    aligned_X, aligned_Y = hirschberg(X[v_start : v_end + 1], Y[w_start : w_end + 1])
    print(aligned_X)
    print(aligned_Y)