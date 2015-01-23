##' Profile Bayesian inversion of PROSPECT
##'
##' Initialize inversion with input parameters as defined in quicktest.sh 
##' Input arguments are the same as those used by runpi.R, which are:
##'   jrsd: Starting SD of Jump distribution; larger values indicate larger jumps
##'   species: Species; must EXACTLY match species tag in spectra filenames
##'   precarg: Whether to use a single precison value for all wavelengths (1) or an individual precision at each wavelength (0)
##'   rearg: Random effects; current options are 'none' and 'leaf'
##'   initarg: How to generate initial conditions; OPTIONS: 'random', 'mle' (R optim function for mean observed spectrum), 'guess' (preset values)
##'   [6] Number of MCMC steps.
##'   [7] Filename tag to identify run results
##'   [8] Sub-folder in 'run_results' for storing output.

#touch run_results/testfolder/empty
#rm run_results/testfolder/*
#Rscript scripts/runpi.R 0.05 Oats 1 leaf random 100 test1 testfolder

# define input arguments
jrsd <- 0.05
species <- 'Oats'
precarg <- 1
if(precarg){
        precision <- 'sp'
} else {
        precision <- 'pwl'
}
rearg <- 'leaf'
initarg <- 'random'
ngibbs <- 10
runid <- 'profile01'
folder <- 'profile'
filename <- sprintf('%s/%s_%g_%s_%s_%s_%s.dat', folder, species, jrsd, precision, rearg, initarg, runid)
dir.create(sprintf('%s', folder), showWarnings = FALSE)

# start profiling
if (!require('proftools', quietly=TRUE)) {
	stop('The proftools package is not installed. Install it!!', call.=FALSE)
}
rprof_out <- sprintf('%s/%s_rprof.txt', folder, runid)
rprof_line <- sprintf('%s/%s_line.txt', folder, runid)
rprof_call <- sprintf('%s/%s_call.txt', folder, runid)
Rprof(filename=rprof_out, interval=0.01, line.profiling=TRUE)

# run
source('inv_bayes.R')
source('specdataproc.R')
smat <- specmatrix(species)
pinvbayes(smat, ngibbs=ngibbs, JumpRSD=jrsd, fname=filename,
          local.store=FALSE, 
          single.precision=precarg, 
          random.effects=rearg,
          inits=initarg,
          ar.step=100)

# stop profiling, generate reports
Rprof(NULL)
sink(rprof_line); summaryRprof(filename=rprof_out, lines='show'); sink()
printProfileCallGraph(readProfileData(rprof_out), file=rprof_call, percent=TRUE)
