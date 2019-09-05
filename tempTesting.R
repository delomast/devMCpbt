### testing script for devMCpbt
# not built as part of package



MCpbt(1, 1, 5, 1)
MCpbt(100, 1000, 7, 10)

rec <- c()
for(i in 1:10000){
	rec <- c(rec,MCpbt(c(1,1,1,1), i)[1])
}
MCpbt(c(1,1,1,1,1), 2)
hist(rec)
summary(rec)

library(MCMCpack)
summary(rdirichlet(10000,c(1,1,1,1))[,1])
rdirichlet(10,c(1,1,1,1,1))
hist(rdirichlet(10000,c(1,1,1,1))[,1])

# storage calculations
iter = 150
burnin = 30
thin = 2
count=0
for(i in 0:(iter - 1)){
	if(i< burnin){next}
	if((i - (burnin -1)) %% thin == 0){
		count = count + 1
	}
}
count
(iter - burnin) %/% thin

# testing clip routine
hist(MCMCclip(100, 500, c(1,1), 15000, 30, 2, 7))

system.time(MCMCclip(100, 500, c(1,1), 15000, 30, 2, 7))
MCMCclip(10, 50, c(1,1), 15000, 30, 2, 7)
