### WTF: TAGA(G)GACT should be equal to TAGA(G)GACT

#chain = 'alpha'
chain = 'beta'
#chain = 'alpha or beta'
#chain = 'alpha and beta'

# alpha clonotypes
WaG_alpha = 'TRAV26-1_AGTCGTTAGTG_TRAJ53'  # alpha
SchH_alpha = 'TRAV12-2_TGAACTGGGATA_TRAJ33' # alpha
FrK_alpha = 'TRAV6_TGTGCCGTACGGGATA_TRAJ12' # alpha

# beta clonotypes
WaG_beta = 'TRBV20-1_CTAGATCGAGACTAGCGGGCTCCT_TRBJ2-7'  # beta
SchH_beta = 'TRBV20-1_AGAGATTTAACTAGCGGGAGCTTCCCCTAC_TRBJ2-7' # beta
FrK_beta = 'TRBV4-1_CCAAGCCCTCGGGGGGGATAAT_TRBJ1-6'  # beta

if chain == 'alpha':
    WaG_clonotype = WaG_alpha
    SchH_clonotype = SchH_alpha
    FrK_clonotype = FrK_alpha
elif chain == 'beta':
    WaG_clonotype = WaG_beta
    SchH_clonotype = SchH_beta
    FrK_clonotype = FrK_beta
elif chain == 'alpha or beta':
    WaG_clonotype = WaG_alpha + ' or ' + WaG_beta
    SchH_clonotype = SchH_alpha + ' or ' + SchH_beta
    FrK_clonotype = FrK_alpha + ' or ' + FrK_beta
elif chain == 'alpha and beta':
    WaG_clonotype = WaG_alpha + ' and ' + WaG_beta
    SchH_clonotype = SchH_alpha + ' and ' + SchH_beta
    FrK_clonotype = FrK_alpha + ' and ' + FrK_beta
else:
    assert(False)



WaG_SE  = [ ('WaG', 'SE',   'sunflower yellow', WaG_clonotype) ] # orig: 'TRBV20-1_CTAGA(G)GACTAGCGGGCTCCT_TRBJ2-1'
WaG_SL  = [ ('WaG', 'SL',   'goldenrod',        WaG_clonotype) ] #'TRAV26-1_AGTCGTTAGTG_TRAJ53 or TRBV20-1_CTAGA(G)GACTAGCGGGCTCCT_TRBJ2-1'), #NO good clones#
WaG_BE  = [ ('WaG', 'BE',   'peach',            WaG_clonotype) ] #'TRAV10_GTGAGGGACCCCATCCTCTG_TRAJ28 or TRBV6-5_GTTACGGGGACAAGTACGA_TRBJ2-7'), #NO good clones#
WaG_BL  = [ ('WaG', 'BL',   'bright blue',      WaG_clonotype) ] #'TRAV35_GGCAGGGGTCAAGTGC_TRAJ3 or TRBV6-5_GTTACGGGGACAAGTACGA_TRBJ2-7'), #NO good clones
WaG_BLC = [ ('WaG', 'BLC',  'dusty rose',       WaG_clonotype) ] # 11 cells
WaG_good = WaG_SE + WaG_BLC
WaG_bad  = WaG_SL + WaG_BE + WaG_BL + WaG_BLC
WaG_all  = WaG_good + WaG_bad

SchH_SE = [ ('SchH', 'SE', 'mustard yellow',  SchH_clonotype) ] #'TRAV12-2_TGAACTGGGATA_TRAJ33 or TRBV20-1_AGAGATTTAACTAGCGGGAGCTTCCCCTAC_TRBJ2-7'),
SchH_SL = [ ('SchH', 'SL', 'dark orange',     SchH_clonotype) ]
SchH_BE = [ ('SchH', 'BE', 'orchid',          SchH_clonotype) ]  # old hypo: 'TRAV13-1_AGCAACCTTTCCTAAC_TRAJ20'
SchH_BL = [ ('SchH', 'BL', 'purple',          SchH_clonotype) ] #NO good clones#
SchH_good = SchH_SE + SchH_SL
SchH_bad  = SchH_BL + SchH_BE
SchH_all  = SchH_good + SchH_bad

FrK_SE = [ ('FrK', 'SE', 'chartreuse', FrK_clonotype) ]
FrK_SL = [ ('FrK', 'SL', 'green',      FrK_clonotype) ] #'TRBV4-1_CCAAGCCCTCGGGGGGGATAAT_TRBJ1-6'),
FrK_BE = [ ('FrK', 'BE', 'coral',      FrK_clonotype) ] #NO good clones#
FrK_BL = [ ('FrK', 'BL', 'crimson',    FrK_clonotype) ] # TRAV6_TGTGCCGTACGGGATA_TRAJ12 or TRBV4-1_CCAAGCCCTCGGGGGGGATAAT_TRBJ1-6'
FrK_good = FrK_SE + FrK_SL + FrK_BL
FrK_bad = FrK_BE
FrK_all = FrK_good + FrK_bad

All_good = WaG_good + SchH_good + FrK_good
All = WaG_all + SchH_all + FrK_all

All_nonpreselected = (WaG_SE + WaG_SL + WaG_BE + WaG_BL) + (FrK_SE + FrK_SL) + (SchH_SE + SchH_SL)
All_preselected = (WaG_BLC) + (FrK_BE + FrK_BL) + (SchH_BE + SchH_BL)

NoClone = 'NoClone'
Tracer_12353 = [ ('tracer', '12353', 'chartreuse', NoClone) ]
Tracer_12392 = [ ('tracer', '12392', 'green', NoClone) ]
Tracer_15123 = [ ('tracer', '15123', 'coral', NoClone) ]
Tracer_15126 = [ ('tracer', '15126', 'crimson', NoClone) ]

fluidigm = [ ('fluidigm', '', 'blue',    NoClone) ]
batch1   = [ ('batch1', '', 'yellow',    NoClone) ]
other_samples = fluidigm + batch1






plate1_LaM_clonotype = 'TRBV20-1_GTGCTCCCCGACTAGCGGGAAGCTTACCCTAC_TRBJ2-7'
plate2_SchU_clonotype = NoClone
plate5_RuW_clonotype = NoClone #'TRBV20-1_GTGCTCCCCGACTAGCGGGAAGCTTACCCTAC_TRBJ2-7'

plate1_LaM  = [ ('Plate1_LaM', '',   'sunflower yellow', plate1_LaM_clonotype) ]
plate2_SchU = [ ('Plate2_SchU', '',  'chartreuse', plate2_SchU_clonotype) ]
plate5_RuW  = [ ('Plate5_RuW', '',   'coral', plate5_RuW_clonotype) ]


Batch5 = plate1_LaM + plate2_SchU + plate5_RuW











