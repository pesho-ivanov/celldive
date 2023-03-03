def get_stats(a_cells, b_cells, a_chains, b_chains, color_f):
    cells_perc  = '{:.3}%'.format(100.0*a_cells/(a_cells + b_cells))
    chains_perc = '{:.3}%'.format(100.0*a_chains/(a_chains + b_chains))
    return '{:2} cells ({}) / {:2} chains ({})  // % of all a&b cells/chains'.format(
        a_cells, color_f(cells_perc), a_chains, color_f(chains_perc))

def count_all(sample):
    print('')
    print('')
    tcrs = sample.tcrs
    alpha_cells, alpha_chains = 0, 0
    beta_cells, beta_chains = 0, 0
    for cell in tcrs.cells():
        alphas, betas = tcrs.cell2tcrs(cell)
        if len(alphas) > 0:
            alpha_cells += 1
            alpha_chains += len(alphas)
        if len(betas) > 0:
            beta_cells += 1
            beta_chains += len(alphas)
    cells = len(tcrs.cells())
    print(sample)
    print('  {}'.format(get_stats(alpha_cells, beta_cells, alpha_chains, beta_chains, green)))
    #print('       alphas: {:2} ({:.3}%) cells / {:2} chains'.format(alpha_cells, 100.0*alpha_cells/cells, alpha_chains))
    #print('        betas: {:2} ({:.3}%) cells / {:2} chains'.format(beta_cells, 100.0*beta_cells/cells, beta_chains))
    #print('a/(a+b) ratio: {:.3}                 {:.3}'.format(
    #        alpha_cells/(alpha_cells+beta_cells), alpha_chains/(alpha_chains+beta_chains)))

def count(sample, which_chain, what_to_count, mute=False):
    hypo = sample.clones.best_hypo()
    tcrs = sample.tcrs
    cl_alpha, noncl_alpha, cl_alpha_sum, noncl_alpha_sum = 0, 0, 0, 0
    cl_beta, noncl_beta, cl_beta_sum, noncl_beta_sum = 0, 0, 0, 0
    for cell in tcrs.cells():
        alphas, betas = tcrs.cell2tcrs(cell)
        if len(alphas) > 0:
            if cell in hypo.clonal_cells:
                cl_alpha += 1
                cl_alpha_sum += len(alphas)
            if cell in hypo.nonclonal_cells:
                noncl_alpha += 1
                noncl_alpha_sum += len(alphas)
        if len(betas) > 0:
            if cell in hypo.clonal_cells:
                cl_beta += 1
                cl_beta_sum += len(betas)
            if cell in hypo.nonclonal_cells:
                noncl_beta += 1
                noncl_beta_sum += len(betas)
        #print('{} -> {}, {}'.format(cell, alphas, betas))
    if not mute:
        cells = len(tcrs.cells())
        print(sample)
        print('      clonal alphas: {}'.format(get_stats(cl_alpha, cl_beta, cl_alpha_sum, cl_beta_sum, red)))
        print('   nonclonal alphas: {}'.format(get_stats(noncl_alpha, noncl_beta, noncl_alpha_sum, noncl_beta_sum, green)))
        
        #print('   nonclonal alphas: {:2} ({:.3}%) cells / {:2} chains'.format(noncl_alpha, 100.0*noncl_alpha/cells, noncl_alpha_sum))
        #print('cl/(cl+noncl) ratio: {:.3}                   {:.3}'.format(
        #    cl_alpha/(cl_alpha+noncl_alpha), cl_alpha_sum/(cl_alpha_sum+noncl_alpha_sum)))
    #scipy.stats.binom_test(cl_alpha, cl_alpha+noncl_alpha, 1.0/2, alternative='greater')
    if what_to_count == 'cells':
        return cl_alpha, noncl_alpha
    elif what_to_count == 'chains':
        return cl_alpha_sum, noncl_alpha_sum
    else:
        assert False

def test_decrease(sample1, sample2, chain):
    print('')
    print('chain:', chain)
    p_vals = []
    for what_to_measure, mute in [('cells', False), ('chains', True)]:
        c1, n1 = count(sample1, chain, 'cells', mute)
        c2, n2 = count(sample2, chain, 'cells', mute)
        p_vals.append( scipy.stats.binom_test(c2, c2+n2, 1.0*c1/(c1+n1), alternative='less') )
        p_vals.append( scipy.stats.binom_test(c1, c1+n1, 1.0*c2/(c2+n2), alternative='greater') ) 
    print('p_vals: [{}, {}]'.format(pval2colstr(min(p_vals)), pval2colstr(max(p_vals))))
    
#test_decrease(run.B['SchH_SE'], run.B['SchH_SL'])
#test_decrease(run.B['FrK_SE'], run.B['FrK_SL'], chain='alpha')
#test_decrease(run.B['FrK_SE'], run.B['FrK_SL'], chain='beta')
#test_decrease(run.B['FrK_BE'], run.B['FrK_BL'])
#test_decrease(run.B['WaG_SE'], run.B['WaG_SL'])
#count_all(run.B['fluidigm'])
