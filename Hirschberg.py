def hirschberg(X, Y):
    """
    Hirschberg algorithm for global alignment.
    :param X: First sequence (string) 
    :param Y: Second (longer) sequence (string) 
    :return: Tuple of two aligned sequences with gaps (X_aligned, Y_aligned)
    """
    def NW_score(X, Y):
        """ 
        Needleman-Wunsch score computation using linear space.
        Returns the last row of the scoring matrix.
        """
        m, n = len(X), len(Y)
        prev = [0] * (n + 1)
        curr = [0] * (n + 1)

        # Initialize first row
        for j in range(1, n + 1):
            prev[j] = prev[j - 1] - 1  # Gap penalty

        for i in range(1, m + 1):
            curr[0] = prev[0] - 1  # Gap penalty
            for j in range(1, n + 1):
                match = prev[j - 1] + (1 if X[i - 1] == Y[j - 1] else -1)  # Match/mismatch
                delete = prev[j] - 1  # Deletion
                insert = curr[j - 1] - 1  # Insertion
                curr[j] = max(match, delete, insert)
            prev, curr = curr, prev

        return prev

    def backtrace(X, Y):
        """Recursive backtrace using the Hirschberg split."""
        if len(X) == 0:
            return "-" * len(Y), Y
        elif len(Y) == 0:
            return X, "-" * len(X)
        elif len(X) == 1 or len(Y) == 1:
            return needleman_wunsch(X, Y)    #base case for needleman-wunsch

        # Split X in half
        mid = len(X) // 2
        score_L = NW_score(X[:mid], Y)  #prefix
        score_R = NW_score(X[mid:][::-1], Y[::-1])[::-1]

        # Find optimal split point
        split = max(range(len(Y) + 1), key=lambda j: score_L[j] + score_R[j])

        # Recursively align the two halves
        X_left, Y_left = backtrace(X[:mid], Y[:split])
        X_right, Y_right = backtrace(X[mid:], Y[split:])
        return X_left + X_right, Y_left + Y_right

    return backtrace(X, Y)

def needleman_wunsch(X, Y):
    """ Basic Needleman-Wunsch for base cases. """
    m, n = len(X), len(Y)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        dp[i][0] = dp[i - 1][0] - 1  # Gap penalty
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j - 1] - 1  # Gap penalty

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i - 1][j - 1] + (1 if X[i - 1] == Y[j - 1] else -1)  # Match/mismatch
            delete = dp[i - 1][j] - 1  # Deletion
            insert = dp[i][j - 1] - 1  # Insertion
            dp[i][j] = max(match, delete, insert)

    # Backtrace
    i, j = m, n
    X_aligned, Y_aligned = "", ""
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + (1 if X[i - 1] == Y[j - 1] else -1):
            X_aligned = X[i - 1] + X_aligned
            Y_aligned = Y[j - 1] + Y_aligned
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] - 1:
            X_aligned = X[i - 1] + X_aligned
            Y_aligned = "-" + Y_aligned
            i -= 1
        else:
            X_aligned = "-" + X_aligned
            Y_aligned = Y[j - 1] + Y_aligned
            j -= 1

    return X_aligned, Y_aligned

# Example usage:
# X = "TGC"
# Y = "ACTGCAA"
# aligned_X, aligned_Y = hirschberg(X, Y)
# print("Aligned Sequences:")
# print(aligned_X)
# print(aligned_Y)
