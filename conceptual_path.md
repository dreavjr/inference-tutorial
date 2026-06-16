The sequence of demonstrations is done, the code is refactored. Now we must revise all text on the pages to make a coherent whole. This is the conceptual parcours:


# Step 1: Classical statistical test

Classical inference doesn't accept parameters as random variables because that would clash with the frequentist notion of probability. To bypass that blocking point, we define the "experiment" as the source of variation. Each experiment is associated to a random sample, thus the importance of sampling distributions and test statistcs. Instead of conceptualizing the uncertainty of a parameter, we will instead conceptualize significance: how often our experiment leads to the right decision (test) or to the right range of possibilities (confidence interval).

We cannot control directly for power, but, under the assumptions of the statistical model, we can compute how much power our current sample size gives us, or vice versa.

Concepts:
- Experimental sample, sampling distribution, test statistic
- Tail areas and p-values
- Statistical tests and confidence intervals

Take-homes:
- Classical inference is easy to apply (at least for traditional problems such as comparing population means) but hard to interpret. Classical tests have lots of assumptions. Significance and power apply across multiple experiments, not one particular experiment
- The bread and butter of classical inference for parameters is maximum-likelihood to locate a parameter and use sampling probabilites and tail-area analysis for determining uncertainty

# Step 2: What we want vs. can control

The ideas of power and significance appear throughout the literature with different names: sensitivity and specificity, true positive rate and true negative rate, etc. They are formal properties of the test but do not correspond to what we will observe in the field, because they disconsider base rates. PPV and NPV correspond much more to what we will see in the field ("the test tagged this as positive/negative: how often can I expect it to be right?") while power and significance answer "the sample is positive/negative: how often can I expect the test to tag it right" (which is not as useful, because the true labels are, of course, unknown: if I had them I wouldn't need a test!).

Take-homes:
- Test-bound knobs are easy (change the test)
- Reality-bound knobs are hard (changing the reality is much harder)
- Sensitivity and specificity don't bound Type I and II errors, not even asymptotically. However, if I make assumptions about the base rate, then I can find useful bounds.

# Step 3: A difference in style?

Bayesian probabilities are not tied to frequencies: they represent what you must conclude when you are reasoning under uncertainty. A probability distribution represent how much you can locate the value of a random variable given all that you know. This small change changes everything, because now parameters don't have to be fixed: their distribution represent what you know about them, not frequencies that depend on repetition.

Concepts:
- Bayesian probabilities break the asymmetry between parameters and data: everything is a random variable with a distribution
- Inference by Bayes' rule: posterior is proportional to prior x likelihood. The posterior represents how much I expect the parameter's value to be in certain regions of the parameter space
- Priors don't have to be informative

Myths to dismiss:
- Bayesian probability is all about priors => you can often use non-informative or little-informative priors (although that is not necessarily a good idea for most problems)
- I can use Bayesian inference to mimick hard-threshold classical tests and still reap the benefits => if you don't add informative priors ("let the data speak for itself" zealousness) and force sharp decisions instead of exploit the uncertainty, Bayesian inference can give exactly the same results as classical inference

Take-homes:
- This is not the full Bayesian treatment: (1) I'm using maximum a posteriori with an uninformative prior = same as maximum likelihood. Credible intervals force a binarization of the decision that goes against the Bayesian ethos.

# Step 4: Bayesian inference's pipeline: prior -> lihelihood -> posterior -> da capo

Full Bayeasian treatment: representation of uncertainty, marginalization of nuisance variables (the variance, in this case). Conclusion is not a single value or hard interval: the posterior _is_ the conclusion.

Concepts:
- Prior, likelihood, posterior

Take-homes:
- Likelihoods are not normalized, nor the product prior x likelihood. Normalizing is _the_ main challenge of Bayesian inference
- Priors should reflect domain knowledge and be lax enough to allow the model to "change opinions", but making them _too_ lax leads to missed opportunities

# Step 5: Everyone is entitled to their opinion (but they'll lose money if they disagree with me)


=> Bibliography
- Accessible resources
- Practical materials and tools
- Academic material



