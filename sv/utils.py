import numpy as np

def editdist(refdist, readdist):
    penalty_delete = 50
    penalty_insert = 50
    # penalty_substitute = 500
    # match_thresh = 30

    reflen = len(refdist)
    readlen = len(readdist)
    matrix = np.zeros((reflen+1, readlen+1))

    refdist_orig = np.copy(refdist)
    readdist_orig = np.copy(readdist)

    for i in range(1, reflen):
        refdist[i] = refdist[i-1] + refdist[i]
    print('refdist: %s' % str(refdist))
    for i in range(1, readlen):
        readdist[i] = readdist[i-1] + readdist[i]
    print('readdist: %s' % str(readdist))

    for i in range(1, reflen+1):
        matrix[i,0] = penalty_delete * i

    for i in range(1, readlen+1):
        matrix[0,i] = penalty_insert * i

    for i in range(1, reflen+1):
        c1 = refdist[i-1]
        for j in range(1, readlen+1):
            c2 = readdist[j-1]
            deltac = abs(c2-c1)

            '''
            if deltac < match_thresh:
                matrix[i,j] = matrix[i-1,j-1] + deltac
            else:
            '''
            if True:
                cost_match = matrix[i-1,j-1] + deltac
                cost_delete = matrix[i-1,j] + penalty_delete
                cost_insert = matrix[i,j-1] + penalty_insert
                # cost_substitute = matrix[i-1,j-1] + penalty_substitute

                cost_minimum = cost_delete
                if cost_insert < cost_minimum:
                    cost_minimum = cost_insert
                # if cost_substitute < cost_minimum:
                    # cost_minimum = cost_substitute
                if cost_match < cost_minimum:
                    cost_minimum = cost_match

                matrix[i,j] = cost_minimum

    print('dist matrix:\n%s' % (str(matrix)))

    # trace backwards
    trace = []
    i = reflen
    j = readlen
    while i > 0 and j > 0:
        # test the 3 cases, pick lowest cost path backwards
        cost_match = matrix[i-1,j-1]
        cost_delete =  matrix[i-1,j]
        cost_insert = matrix[i,j-1]

        if cost_match < cost_delete and cost_match < cost_insert:
            i -= 1
            j -= 1
            trace.insert(0, 'match')
        elif cost_delete < cost_insert:
            i -= 1
            trace.insert(0, 'delete')
        else:
            j -= 1
            trace.insert(0, 'insert')

    print('trace:\n')
    for i in range(len(trace)):
        print(' %s' % trace[i])

    return matrix[reflen,readlen]

if __name__ == "__main__":
    ref = [200, 225, 185, 202, 178, 198, 232, 217, 206]
    read1 = [185, 202, 178, 198, 232]

    print('comparing read 1')
    dist = editdist(ref[2:2+5], read1)
    print('dist: %f' % dist)

    print('comparing read 2')
    read2 = [189, 200, 187, 202, 230]
    dist = editdist(ref[2:2+5], read2)
    print('dist: %f' % dist)

    print('comparing read 3')
    read3 = [189, 200, 187, 432]
    dist = editdist(ref[2:2+5], read3)
    print('dist: %f' % dist)

    print('comparing read 4')
    read3 = [189, 200, 90, 97, 202, 230]
    dist = editdist(ref[2:2+5], read3)
    print('dist: %f' % dist)

