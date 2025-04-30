
    cov_matrices = [
        inv(visp_f_reduced),
        inv(visb_f_reduced),
        inv(visb_f_reduced + visp_f_reduced)
    ]

    labels = ['powerspectr