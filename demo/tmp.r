
# Generate data
## Set up SIR model
sckm = list()
sckm$s = 3 # species (S,I,R)
sckm$r = 2 # reactions (S->I, I->R)
#                   S -> I    I -> R
sckm$Pre  = rbind( c(1,1,0), c(0,1,0))
sckm$Post = rbind( c(0,2,0), c(0,0,1))
sckm$stoich = t(sckm$Post-sckm$Pre)
sckm$X = c(16000,100,0)
N = sum(sckm$X)
sckm$theta = c(0.5,0.25)
sckm$lmult = log(c(1/N,1))

## Simulate data
set.seed(2)
n = 50

### True states and transitions
tl = tau.leap(sckm, n)

### Sample transitions
p = c(0.5,0.5) # Sample probabilities for S->I and I->R respectively
y = t(cbind(rbinom(n, tl$nr[,1], p[1]), rbinom(n, tl$nr[,2], p[2])))

# Perform inference
cat("Running sequential inference...\n")
prior = tlpl.prior(sckm$X, 1e1, 1e1, sckm$theta*2, 2, sckm$r)

it = sample(1e6,1)
set.seed(it)
z = tlpl(list(y=y, tau=1), sckm, prior=prior, n.particles=1e3, engine="R", verbose=1)
R.state = .Random.seed


set.seed(it)
z = tlpl(list(y=y, tau=1), sckm, prior=prior, n.particles=1e3, engine="C", verbose=1)
C.state = .Random.seed

if (!isTRUE(all.equal(R.state, C.state))) print("Ending states are not the same.\n")

