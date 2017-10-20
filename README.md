**PROBLEM STATEMENT**

[The Gene Expression Omnibus (GEO) data series GSE4115 contains data
from 192 human subjects, each with 22,283 profiled genes. Each subject
can have one of three disease states: cancer, no cancer, or suspected
cancer](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4115). Your
task is to build a classifier for cancer vs. no cancer by using HDLSS
techniques (such as elastic net).

**NOTES**\
For the purpose of classification into 2 classes a.k.a cancer or no
cancer, the subjects numbered 188 to 192 are ignored in building the
model. Upon analysing the data, it is found that in each of the last 67
features, more than 75% of the information is marked ‘not available’.
Thus, as a pre-processing step, those features are not considered in
building the model. For the first three models discussed below, the
penalty is defined as (lambda\*((1 − α)\*|β|~2 + α\*|β|~1)) where
alpha=1 corresponds to the lasso penalty, alpha=0 to the ridge penalty
and the elastic net corresponding to 0 ≤ α ≤ 1. The number of lambda
values considered for each of the above models is 100, the default value
in *glmnet()* package.

Cross-validation technique is used in building the models and one of the
main reasons for using cross-validation instead of using the
conventional validation (e.g. partitioning the data set into two sets of
70% for training and 30% for test) is that there is not enough data
available to partition it into separate training and test sets without
losing significant modelling or testing capability. Thus, in this case,
a fair way to properly estimate model prediction performance is to use
cross-validation as a powerful general technique. In summary,
cross-validation combines (averages) measures of fit (prediction error)
to derive a more accurate estimate of model prediction performance.

**refer to Readme.pdf for detailed explanation with examples and images**


