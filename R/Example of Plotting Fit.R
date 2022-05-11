# Plot: SR Figures =============================================================

pdf(file=file.path(dir.figs, "Ricker Parameters.pdf"), height=5, width=6)
par(mfrow=c(1,3), oma=c(1,1,3,1), mar=c(4,1,0,1))

# Ricker Parameters
plotPost(pars$ricker_alpha, xlab="Alpha")
plotPost(pars$ricker_beta, xlab="Beta")
plotPost(pars$ricker_sigma, xlab="Sigma")
mtext("Ricker Parameters", side=3, outer=TRUE, font=2)
# dev.off()

# SR Fit

par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(4,4,3,1))

quant.rec <- apply(pars$pred_rec, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
ord <- order(spawn)
plot(x=spawn/1e3, y=rec/1e3, type="p", pch=21, bg="blue",
     xlab="Spawning Abundance (thousands)",
     ylab="Recruitment (thousands)",
     ylim=c(0,max(quant.rec,rec))/1e3,
     main="Ricker Model Fit")
grid(col="black")
# Predicted
polygon(x=c(spawn[ord], rev(spawn[ord]))/1e3,
        y=c(quant.rec[1,ord], rev(quant.rec[5,ord]))/1e3,
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x=c(spawn[ord], rev(spawn[ord]))/1e3,
        y=c(quant.rec[2,ord], rev(quant.rec[4,ord]))/1e3,
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x=spawn[ord]/1e3, quant.rec[3,ord]/1e3, col="red", lwd=2)

# Simulate Complete series
trial.spawn <- seq(from=0, to=max(spawn), length.out=100)
n.trial.spawn <- length(trial.spawn)

sim.rec <- matrix(nrow=n.sims, ncol=n.trial.spawn)
i <- 1
for(i in 1:n.trial.spawn) {
  sim.rec[,i] <- trial.spawn[i]*exp(pars$ricker_alpha - pars$ricker_beta*trial.spawn[i])
}
quant.sim.rec <- apply(sim.rec, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

# Plot Complete series
plot(x=spawn/1e3, y=rec/1e3, type="p", pch=21, bg="blue",
     xlab="Spawning Abundance (thousands)",
     ylab="Recruitment (thousands)",
     ylim=c(0,max(quant.sim.rec,rec))/1e3,
     xlim=c(min(trial.spawn),max(trial.spawn))/1e3,
     main="Ricker Model Fit")
grid(col="black")
# Predicted
polygon(x=c(trial.spawn, rev(trial.spawn))/1e3,
        y=c(quant.sim.rec[1,], rev(quant.sim.rec[5,]))/1e3,
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x=c(trial.spawn, rev(trial.spawn))/1e3,
        y=c(quant.sim.rec[2,], rev(quant.sim.rec[4,]))/1e3,
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x=trial.spawn/1e3, quant.sim.rec[3,]/1e3, col="red", lwd=2)

# Plot SR Prediction
plotPost(pars$curr_fcst_sr, xlab="SR Forecast")

dev.off()