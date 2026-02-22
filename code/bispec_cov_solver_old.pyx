from libc.stdio cimport printf

# the whole point of this file is to calculate C^(-1)C^(-1)C^(-1)B, where C is the covariance matrix and B is the bispectrum vector. we want to do this for all combinations of l_i, i.e. for all XYZ. we want to do this in a way that minimizes the number of times we have to apply C^(-1), since that is the most expensive part. so we will first apply C^(-1) along X, then along Y, then along Z. this way we only have to apply C^(-1) 5*5 = 25 times along X, then 5*5 = 25 times along Y, then 5*5 = 25 times along Z, for a total of 75 applications of C^(-1). if we were to apply C^(-1) directly to B for each XYZ, we would have to apply it 125 times.

cdef inline int index_3d_5(int X, int Y, int Z) noexcept nogil:
    return X * 25 + Y * 5 + Z

cdef inline int index_2d_5(int X, int Y) noexcept nogil:
    return X * 5 + Y

cdef int soe_solver_5(double* matrix, double* vec, double* sol) noexcept nogil:
    cdef int i = 0
    cdef int j = 0
    cdef int k
    cdef int pivot_row
    cdef double maxabs
    cdef double[25] A
    cdef double[5] b
    cdef double[5] y
    cdef double tmp
    cdef double piv


    while i  < 25:
        A[i] = matrix[i]
        i+=1
    i=0
    
    while i < 5:
        b[i] = vec[i]
        i+=1
    i=0


    while i < 5:
        pivot_row = i
        maxabs = A[index_2d_5(i, i)]
        if maxabs < 0.0:
            maxabs = -1*maxabs
        
        # find max abs val in ith collumn below the diagonal
        j=i+1
        while j < 5:

            tmp = A[index_2d_5(j, i)]
            if tmp < 0.0:
                tmp = -tmp
            if tmp > maxabs:
                maxabs = tmp
                pivot_row = j
            j+=1
            # TODO: test for singularity?
        j=0

        # swap current row with pivot row in A and b
        if pivot_row != i:
            # j goes over the elements within a row here
            while j < 5:
                tmp = A[index_2d_5(i, j)]
                A[index_2d_5(i, j)] = A[index_2d_5(pivot_row, j)]
                A[index_2d_5(pivot_row, j)] = tmp
                j+=1
            j=0

            tmp = b[i]
            b[i] = b[pivot_row]
            b[pivot_row] = tmp
        
        # eliminate entries below pivot
        piv = A[index_2d_5(i, i)]
        j = i+1
        while j < 5:
            A[index_2d_5(j, i)] /= piv # L value (A contains L and U at the same time)
            tmp = A[index_2d_5(j, i)]
            
            # subtract the rows
            k = i+1
            while k < 5:
                A[index_2d_5(j, k)] -= tmp * A[index_2d_5(i, k)]
                k += 1
            k = 0
            j+=1
        j=0
        i+=1
    i=0

    # Forward substitution, i.e. solve L y = b
    while i < 5:
        tmp = b[i]
        while j < i:
            tmp -= A[index_2d_5(i, j)] * y[j]
            j+=1
        j=0
        y[i] = tmp

        i+=1
    i=0

    i=4
    while i > -1:
        tmp = y[i]
        
        j=i+1
        while j < 5:
            tmp -= A[index_2d_5(i, j)] * sol[j]

            j+=1
        j=0

        piv = A[index_2d_5(i, i)]
        sol[i] =  tmp / piv
        
        i-=1

    return 0

cdef int solve_1(double* CXp, int Yp, int Zp, double* bispec_vec, double* out5) noexcept nogil:
    # here "bispec vec" refers to all the values of (d/d alpha) B^(XYZ)_(l_1l_2l_3) for fixed alpha and l_i, i.e. bispec_vec[index_3d(X,Y,Z)]
    cdef double[5] vec
    cdef int Xp = 0
    while Xp < 5:
        vec[Xp] = bispec_vec[index_3d_5(Xp, Yp, Zp)]
        Xp += 1

    soe_solver_5(CXp, &vec[0], out5) # so this will be a vector. for given XYZ this function is just a function of Yp and Zp
    return 0

cdef int solve_2(double* CXp, double* CYp, int Zp, double* bispec_vec, double* out25) noexcept nogil:
    cdef double[25] vecs
    cdef double[5] b
    cdef double[5] x
    cdef int Xp = 0
    cdef int Yp = 0

    # build vecs[Xp,Yp], (C1^{-1}_{l1})^XXp applied on Xp
    Yp=0
    while Yp < 5:
        solve_1(CXp, Yp, Zp, bispec_vec, &x[0]) # in python language: x = solve_1(...)
        Xp=0
        while Xp < 5:
            vecs[index_2d_5(Xp, Yp)] = x[Xp]
            Xp += 1
        Yp += 1
        
    # for each Xp, solve along Y: Cy * y = vecs[Xp,:]
    Xp=0
    while Xp < 5:
        Yp=0
        while Yp < 5:
            b[Yp] = vecs[index_2d_5(Xp, Yp)]
            Yp += 1
        soe_solver_5(CYp, &b[0], &x[0])
        Yp=0
        while Yp < 5:
            out25[index_2d_5(Xp, Yp)] = x[Yp]
            Yp += 1
        Xp += 1

    return 0

cdef int solve_3(double* CXp, double* CYp, double* CZp, double* bispec_vec, double* out125) noexcept nogil:
    cdef double[125] vecss
    cdef double[5] b
    cdef double[5] x
    cdef double[25] x2
    cdef int Xp = 0
    cdef int Yp = 0
    cdef int Zp = 0

    # build vecs[Xp,Yp], (C1^{-1}_{l1})^XXp applied on Xp
    Zp = 0
    while Zp < 5:
        solve_2(CXp, CYp, Zp, bispec_vec, &x2[0]) # in python language: x2 = solve_2(...)
        Yp = 0
        while Yp < 5:
            Xp = 0
            while Xp < 5:
                vecss[index_3d_5(Xp, Yp, Zp)] = x2[index_2d_5(Xp, Yp)]
                Xp += 1
            Yp += 1
        Zp += 1

    # for each Xp, Yp, solve along Z: Cz * z = vecs[Zp,:]
    Xp = 0
    while Xp < 5:
        Yp = 0
        while Yp < 5:
            Zp = 0
            while Zp < 5:
                b[Zp] = vecss[index_3d_5(Xp, Yp, Zp)]
                Zp += 1
            soe_solver_5(CZp, &b[0], &x[0])
            Zp = 0
            while Zp < 5:
                out125[index_3d_5(Xp, Yp, Zp)] = x[Zp]
                Zp += 1
            Yp += 1
        Xp += 1

    return 0

def soe_solver_5_wrap(double[::1] matrix, double[::1] vec, double[::1] sol):
    soe_solver_5(&matrix[0], &vec[0], &sol[0])
    return None

def solve_1_wrap(double[::1] CXp, int Yp, int Zp, double[::1] bispec_vec, double[::1] out5):
    solve_1(&CXp[0], Yp, Zp, &bispec_vec[0], &out5[0])
    return None

def solve_2_wrap(double[::1] CXp, double[::1] CYp, int Zp, double[::1] bispec_vec, double[::1] out25):
    solve_2(&CXp[0], &CYp[0], Zp, &bispec_vec[0], &out25[0])
    return None

def solve_3_wrap(double[::1] CXp, double[::1] CYp, double[::1] CZp, double[::1] bispec_vec, double[::1] out125):
    solve_3(&CXp[0], &CYp[0], &CZp[0], &bispec_vec[0], &out125[0])
    return None
