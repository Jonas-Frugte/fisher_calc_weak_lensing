from libc.stdio cimport printf

# Generalized for variable number of tracers (n_tracers)
cdef int index_3d(int X, int Y, int Z, int n_tracers) noexcept nogil:
    return X * n_tracers * n_tracers + Y * n_tracers + Z

cdef int index_2d(int X, int Y, int n_tracers) noexcept nogil:
    return X * n_tracers + Y

cdef int soe_solver(
    double* matrix, # n_tracer^2
    double* vec, # n_tracer
    double* sol, # n_tracer
    double* A, # n_tracer^2
    double* b, # n_tracer
    double* y, # n_tracer
    int n_tracers
    ) noexcept nogil:
    cdef int i = 0
    cdef int j = 0
    cdef int k
    cdef int pivot_row
    cdef double maxabs
    cdef double tmp
    cdef double piv


    while i  < n_tracers**2:
        A[i] = matrix[i]
        i+=1
    i=0
    
    while i < n_tracers:
        b[i] = vec[i]
        i+=1
    i=0


    while i < n_tracers:
        pivot_row = i
        maxabs = A[index_2d(i, i, n_tracers)]
        if maxabs < 0.0:
            maxabs = -1*maxabs
        
        # find max abs val in ith collumn below the diagonal
        j=i+1
        while j < n_tracers:

            tmp = A[index_2d(j, i, n_tracers)]
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
            while j < n_tracers:
                tmp = A[index_2d(i, j, n_tracers)]
                A[index_2d(i, j, n_tracers)] = A[index_2d(pivot_row, j, n_tracers)]
                A[index_2d(pivot_row, j, n_tracers)] = tmp
                j+=1
            j=0

            tmp = b[i]
            b[i] = b[pivot_row]
            b[pivot_row] = tmp
        
        # eliminate entries below pivot
        piv = A[index_2d(i, i, n_tracers)]
        j = i+1
        while j < n_tracers:
            A[index_2d(j, i, n_tracers)] /= piv # L value (A contains L and U at the same time)
            tmp = A[index_2d(j, i, n_tracers)]
            
            # subtract the rows
            k = i+1
            while k < n_tracers:
                A[index_2d(j, k, n_tracers)] -= tmp * A[index_2d(i, k, n_tracers)]
                k += 1
            k = 0
            j+=1
        j=0
        i+=1
    i=0

    # Forward substitution, i.e. solve L y = b
    while i < n_tracers:
        tmp = b[i]
        while j < i:
            tmp -= A[index_2d(i, j, n_tracers)] * y[j]
            j+=1
        j=0
        y[i] = tmp

        i+=1
    i=0

    i = n_tracers - 1
    while i > -1:
        tmp = y[i]
        
        j=i+1
        while j < n_tracers:
            tmp -= A[index_2d(i, j, n_tracers)] * sol[j]

            j+=1
        j=0

        piv = A[index_2d(i, i, n_tracers)]
        sol[i] =  tmp / piv
        
        i-=1

    return 0

cdef int solve_1(
    double* CXp, # n_tracer^2
    int Yp, 
    int Zp, 
    double* bispec_vec, # n_tracer^3
    double* outn, # n_tracer
    double* vec, # n_tracer 
    double* A_solve, # n_tracer^2
    double* b_solve, # n_tracer
    double* y_solve, # n_tracer
    int n_tracers
    ) noexcept nogil:
    # here "bispec vec" refers to all the values of (d/d alpha) B^(XYZ)_(l_1l_2l_3) for fixed alpha and l_i, i.e. bispec_vec[index_3d(X,Y,Z)]
    cdef int Xp = 0
    while Xp < n_tracers:
        vec[Xp] = bispec_vec[index_3d(Xp, Yp, Zp, n_tracers)]
        Xp += 1

    soe_solver(CXp, vec, outn, A_solve, b_solve, y_solve, n_tracers) # so this will be a vector. for given XYZ this function is just a function of Yp and Zp
    return 0

cdef int solve_2(
    double* CXp, # n_tracer^2
    double* CYp, # n_tracer^2
    int Zp, 
    double* bispec_vec, # n_tracer^3
    double* outn2, # n_tracer^2
    double* vecs, # n_tracer^2
    double* b, # n_tracer
    double* x, # n_tracer
    double* vec_solve1, # n_tracer
    double* A_solve, # n_tracer^2
    double* b_solve, # n_tracer
    double* y_solve, # n_tracer
    int n_tracers
    ) noexcept nogil:

    cdef int Xp = 0
    cdef int Yp = 0

    # build vecs[Xp,Yp], (C1^{-1}_{l1})^XXp applied on Xp
    Yp=0
    while Yp < n_tracers:
        solve_1(CXp, Yp, Zp, bispec_vec, x, vec_solve1, A_solve, b_solve, y_solve, n_tracers) # in python language: x = solve_1(...)
        Xp=0
        while Xp < n_tracers:
            vecs[index_2d(Xp, Yp, n_tracers)] = x[Xp]
            Xp += 1
        Yp += 1
        
    # for each Xp, solve along Y: Cy * y = vecs[Xp,:]
    Xp=0
    while Xp < n_tracers:
        Yp=0
        while Yp < n_tracers:
            b[Yp] = vecs[index_2d(Xp, Yp, n_tracers)]
            Yp += 1
        soe_solver(CYp, b, x, A_solve, b_solve, y_solve, n_tracers)
        Yp=0
        while Yp < n_tracers:
            outn2[index_2d(Xp, Yp, n_tracers)] = x[Yp]
            Yp += 1
        Xp += 1

    return 0

cdef int solve_3(
    double* CXp, # n_tracer^2
    double* CYp, # n_tracer^2
    double* CZp, # n_tracer^2
    double* bispec_vec, # n_tracer^3
    double* outn3, # n_tracer^3
    double* vecss, # n_tracer^3
    double* b, # n_tracer
    double* x, # n_tracer
    double* x2, # n_tracer^2
    double* vecs_solve2, # n_tracer^2
    double* b_solve2, # n_tracer
    double* x_solve2, # n_tracer
    double* vec_solve1, # n_tracer
    double* A_solve, # n_tracer^2
    double* b_solve, # n_tracer
    double* y_solve, # n_tracer
    int n_tracers
    ) noexcept nogil:
    cdef int Xp = 0
    cdef int Yp = 0
    cdef int Zp = 0

    # build vecs[Xp,Yp], (C1^{-1}_{l1})^XXp applied on Xp
    Zp = 0
    while Zp < n_tracers:
        solve_2(CXp, CYp, Zp, bispec_vec, x2, vecs_solve2, b_solve2, x_solve2, vec_solve1, A_solve, b_solve, y_solve, n_tracers) # in python language: x2 = solve_2(...)
        Yp = 0
        while Yp < n_tracers:
            Xp = 0
            while Xp < n_tracers:
                vecss[index_3d(Xp, Yp, Zp, n_tracers)] = x2[index_2d(Xp, Yp, n_tracers)]
                Xp += 1
            Yp += 1
        Zp += 1

    # for each Xp, Yp, solve along Z: Cz * z = vecs[Zp,:]
    Xp = 0
    while Xp < n_tracers:
        Yp = 0
        while Yp < n_tracers:
            Zp = 0
            while Zp < n_tracers:
                b[Zp] = vecss[index_3d(Xp, Yp, Zp, n_tracers)]
                Zp += 1
            soe_solver(CZp, b, x, A_solve, b_solve, y_solve, n_tracers)
            Zp = 0
            while Zp < n_tracers:
                outn3[index_3d(Xp, Yp, Zp, n_tracers)] = x[Zp]
                Zp += 1
            Yp += 1
        Xp += 1

    return 0