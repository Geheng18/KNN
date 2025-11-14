#library(Matrix)  # For nearPD function
#library(glmnet)    # For glmnet function (used in KNN.numerical)
#library(CompQuadForm)
#' @importFrom stats coef dnbinom glm lm residuals rnorm var
#' @importFrom Matrix nearPD
#' @importFrom CompQuadForm liu davies
#' @importFrom glmnet glmnet cv.glmnet


#' @title KNN Mixture Chi-Square Overall Test
#'
#' @description
#' Performs an overall test for variance components in a KNN model using a mixture of chi-square distributions.
#'
#' @param phenotype Numeric vector of response variable (phenotype) of length N.
#' @param genotype List of N x p genotype matrices.
#' @param test.index Indices of variance components to test.
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param MINQUE.type Type of MINQUE estimator ('MINQUE0' or 'MINQUE1') (default is 'MINQUE0').
#' @param KNN.type Type of KNN model ('KNN', 'LMM').
#' @param lambda Regularization parameter (default is 0).
#' @param weight Optional weight vector for variance components.
#' @param fixed.effects Optional fixed effects design matrix.
#' @param constrain Logical indicating whether to constrain variance components to be non-negative (default is FALSE).
#'
#' @return A list containing the p-value of the test.
#' @export
KNN.overall.test = function(phenotype, genotype, test.index, input.kernel = NULL,
                            MINQUE.type='MINQUE0', KNN.type='KNN', lambda=0, weight=NULL,
                            fixed.effects=NULL, constrain=FALSE) {
  # Estimate variance components using KNN
  result = KNN(phenotype=phenotype, genotype=genotype, input.kernel=input.kernel,
               MINQUE.type=MINQUE.type, KNN.type=KNN.type, lambda=lambda,
               weight=weight, fixed.effects=fixed.effects, constrain=constrain)
  covariance.matrix = result$covariance.matrix
  theta = result$theta
  precision.matrix = result$precision.matrix
  variance.component.list = result$variance.component.list

  # Compute sigma under H1
  sigma.hypothesis.one = Reduce(f = "+", x = mapply(FUN = "*", variance.component.list,
                                                    covariance.matrix[1,], SIMPLIFY = FALSE))

  # Compute sigma under H2
  sigma.hypothesis.two.list = list()
  test.weight = rep(1, length(test.index))
  index.counter = 1
  for (component.index in test.index) {
    sigma.hypothesis.two = Reduce(f = "+", x = mapply(FUN = "*", variance.component.list,
                                                      covariance.matrix[component.index,], SIMPLIFY = FALSE))
    sigma.hypothesis.two.list[[index.counter]] = sigma.hypothesis.two
    test.weight[index.counter] = covariance.matrix[component.index, component.index]
    index.counter = index.counter + 1
  }

  ratio = sum(theta[test.index] * test.weight)
  sigma.hypothesis.two = Reduce(f = "+", x = mapply(FUN = "*", sigma.hypothesis.two.list,
                                                    test.weight, SIMPLIFY = FALSE))
  ratio = ratio / theta[1]

  # Compute total sigma
  sigma.total = precision.matrix %*% (sigma.hypothesis.two - ratio * sigma.hypothesis.one) %*% precision.matrix

  # Compute eigenvalues and p-value using liu method
  eigenvalues = eigen(sigma.total, symmetric=TRUE, only.values=TRUE)$values

  # Remove very small eigenvalues for numerical stability
  eigenvalues = eigenvalues[eigenvalues > 1e-15]

  p.value = tryCatch({
    davies(q = 0, lambda = eigenvalues)$Qq
  }, error = function(e) {
    # If Liu fails, use Davies
    liu(q = 0, lambda = eigenvalues)
  })

  return(p.value)
}

#' @title KNN Mixture Chi-Square Individual Test
#'
#' @description
#' Performs an individual test for a variance component in a KNN model using a mixture of chi-square distributions.
#'
#' @param phenotype Numeric vector of response variable (phenotype) of length N.
#' @param genotype List of N x p genotype matrices.
#' @param test.index Index of the variance component to test.
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param MINQUE.type Type of MINQUE estimator ('MINQUE0' or 'MINQUE1') (default is 'MINQUE0').
#' @param KNN.type Type of KNN model ('KNN', 'LMM').
#' @param lambda Regularization parameter (default is 0).
#' @param weight Optional weight vector for variance components.
#' @param fixed.effects Optional fixed effects design matrix.
#' @param constrain Logical indicating whether to constrain variance components to be non-negative (default is FALSE).
#'
#' @return A list containing the p-value of the test.
#' @export
KNN.ind.test = function(phenotype, genotype, test.index, input.kernel = NULL,
                        MINQUE.type='MINQUE0', KNN.type='KNN', lambda=0, weight=NULL,
                        fixed.effects=NULL, constrain=FALSE) {
  # Estimate variance components under the full model
  result.full = KNN(phenotype=phenotype, genotype=genotype, input.kernel=input.kernel,
                    MINQUE.type=MINQUE.type, KNN.type=KNN.type, lambda=lambda,
                    weight=weight, fixed.effects=fixed.effects, constrain=constrain)
  covariance.matrix = result.full$covariance.matrix
  theta.full = result.full$theta
  precision.matrix = result.full$precision.matrix
  variance.component.list.full = result.full$variance.component.list

  # Remove the variance component to test to get the reduced model
  variance.component.list.reduced = variance.component.list.full[-test.index]

  # Estimate variance components under the reduced model
  result.reduced = KNN(phenotype=phenotype, genotype=genotype, input.kernel=input.kernel,
                       MINQUE.type=MINQUE.type, KNN.type=KNN.type, lambda=lambda,
                       weight=weight, fixed.effects=fixed.effects, constrain=constrain,
                       variance.component.list=variance.component.list.reduced)
  theta.reduced = result.reduced$theta

  # Compute sigma H1 and sigma H2
  sigma.hypothesis.one = Reduce(f = "+", x = mapply(FUN = "*", variance.component.list.reduced,
                                                    theta.reduced, SIMPLIFY=FALSE))
  sigma.hypothesis.two = Reduce(f = "+", x = mapply(FUN = "*", variance.component.list.full,
                                                    covariance.matrix[test.index,], SIMPLIFY=FALSE))

  # Compute sigma H
  sigma.h1.svd = svd(sigma.hypothesis.one)
  sigma.h1.sqrt = sigma.h1.svd$u %*% diag(sqrt(sigma.h1.svd$d)) %*% t(sigma.h1.svd$v)
  sigma.hypothesis = sigma.h1.sqrt %*% precision.matrix %*% sigma.hypothesis.two %*%
    precision.matrix %*% sigma.h1.sqrt

  # Compute eigenvalues and p-value using liu method
  eigenvalues = eigen(sigma.hypothesis, symmetric=TRUE, only.values=TRUE)$values

  # Remove very small eigenvalues for numerical stability
  eigenvalues = eigenvalues[eigenvalues > 1e-15]

  p.value = tryCatch({
    davies(q=theta.full[test.index], lambda=eigenvalues)$Qq
  }, error = function(e) {
    # If Liu fails, use Davies
    liu(q=theta.full[test.index], lambda=eigenvalues)
  })

  return(p.value)
}


############################################################
# Combined KNN Estimation Function with Integrated CV
############################################################

#' @title KNN Estimation with Integrated Cross-Validation
#'
#' @description
#' Estimates variance components in a Penalized Kernel Neural Network (KNN) model using MINQUE.
#' If lambda is not specified, performs cross-validation to select optimal lambda.
#'
#' @param phenotype Numeric vector of response variable (phenotype) of length N.
#' @param genotype List of N x p genotype matrices.
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param MINQUE.type Type of MINQUE estimator ('MINQUE0' or 'MINQUE1') (default is 'MINQUE0').
#' @param KNN.type Type of KNN model ('KNN', 'LMM')(default is 'KNN').
#' @param lambda Regularization parameter. If NULL, cross-validation will be performed.
#' @param lambda.list List of lambda values to try in cross-validation (used only if lambda is NULL).
#' @param weight Optional weight vector for variance components.
#' @param fixed.effects Optional fixed effects design matrix.
#' @param constrain Logical indicating whether to constrain variance components to be non-negative (default is FALSE).
#' @param cv.criteria Criterion for selecting optimal lambda in CV ('LOSS', 'MSE', 'COR') (default is 'LOSS').
#' @param variance.component.list Optional precomputed list of variance components.
#'
#' @return A list containing estimated variance components (`theta`), fixed effects coefficients (`beta`),
#'         C matrix (`covariance.matrix`), Q matrix (`precision.matrix`), variance component list, and optimal lambda (`lambda.opt`) if CV was performed.
#' @export
KNN = function(phenotype, genotype, input.kernel = NULL, MINQUE.type='MINQUE0', KNN.type='KNN',
               lambda=NULL, lambda.list=NULL, weight=NULL, fixed.effects=NULL, constrain=FALSE,
               cv.criteria='COR', variance.component.list=NULL) {

  # If lambda is not specified, perform cross-validation
  if (is.null(lambda)) {
    # Initialize lambda list if not provided
    if (is.null(lambda.list)) {
      lambda.list = c(0.05, 5, 500, 50000, 5000000)
    }

    # Perform cross-validation
    lambda.opt = perform.cv(phenotype=phenotype, genotype=genotype, input.kernel=input.kernel,
                            MINQUE.type=MINQUE.type, KNN.type=KNN.type,
                            lambda.list=lambda.list, weight=weight, fixed.effects=fixed.effects,
                            constrain=constrain, criteria=cv.criteria)

    # Use the optimal lambda for the final fit
    lambda = lambda.opt
  } else {
    lambda.opt = NULL
  }

  # If variance components are not provided, compute them
  if (is.null(variance.component.list)) {
    variance.component.list = KNN2LMM(genotype=genotype, input.kernel=input.kernel, KNN.type=KNN.type)
  }

  n.variance.component = length(variance.component.list)
  sample.size = length(phenotype)

  # Initialize V matrix based on MINQUE type
  if (MINQUE.type == 'MINQUE0') {
    variance.matrix = diag(1, sample.size, sample.size)
  } else if (MINQUE.type == 'MINQUE1') {
    if (is.null(weight)) {
      weight = rep(1 / n.variance.component, n.variance.component)
    }
    variance.matrix = Reduce(`+`, mapply(`*`, variance.component.list, weight, SIMPLIFY = FALSE))
  }

  # Ensure V matrix is positive definite
  npd.result = nearPD(variance.matrix)
  variance.matrix = as.matrix(npd.result$mat)
  precision.matrix = solve(variance.matrix)

  # Compute R*variance and R*phenotype
  rotated.variance.list = list()
  if (is.null(fixed.effects)) {
    beta = NULL
  } else {
    linear.model = lm(phenotype ~ fixed.effects - 1)
    beta = as.numeric(coef(linear.model))
    phenotype = residuals(linear.model)
  }
  for (component.index in seq.int(n.variance.component)) {
    rotated.variance.list[[component.index]] = precision.matrix %*% variance.component.list[[component.index]]
  }
  rotated.phenotype = precision.matrix %*% phenotype

  # Compute U matrix
  u.matrix = numeric(n.variance.component)
  for (component.index in seq.int(n.variance.component)) {
    u.matrix[component.index] = sum(rotated.phenotype * (variance.component.list[[component.index]] %*% rotated.phenotype))
  }

  # Compute C matrix
  covariance.matrix = matrix(0, n.variance.component, n.variance.component)
  for (row.index in seq.int(n.variance.component)) {
    for (col.index in seq.int(n.variance.component)) {
      covariance.matrix[row.index, col.index] = sum(rotated.variance.list[[row.index]] * rotated.variance.list[[col.index]])
    }
  }

  # Adjust C matrix with penalty
  penalty.matrix = diag(lambda, n.variance.component, n.variance.component)
  penalty.matrix[1, 1] = 0  # No penalty on the first variance component
  covariance.matrix = covariance.matrix + penalty.matrix

  # Ensure C matrix is positive definite
  npd.result = nearPD(covariance.matrix)
  covariance.matrix = as.matrix(npd.result$mat)

  # Estimate variance components
  covariance.matrix = solve(covariance.matrix)
  theta = covariance.matrix %*% u.matrix

  theta = as.numeric(theta)
  if (constrain) {
    theta[theta < 0] = 0
  }

  # Return results including lambda.opt if CV was performed
  result = list(theta=theta, beta=beta, covariance.matrix=covariance.matrix,
                precision.matrix=precision.matrix, variance.component.list=variance.component.list)

  if (!is.null(lambda.opt)) {
    result$lambda.opt = lambda.opt
  }

  return(result)
}

############################################################
# Internal Cross-Validation Function
############################################################

#' @title Internal Cross-Validation Function
#'
#' @description
#' Internal function to perform cross-validation for lambda selection.
#' This function is called internally by KNN when lambda is not specified.
#'
#' @param phenotype Numeric vector of response variable (phenotype) of length N.
#' @param genotype List of N x p genotype matrices.
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param MINQUE.type Type of MINQUE estimator ('MINQUE0' or 'MINQUE1').
#' @param KNN.type Type of KNN model ('KNN', 'LMM').
#' @param lambda.list List of lambda values to try.
#' @param weight Optional weight vector for variance components.
#' @param fixed.effects Optional fixed effects design matrix.
#' @param constrain Logical indicating whether to constrain variance components to be non-negative.
#' @param criteria Criterion for selecting optimal lambda ('MSE', 'COR').
#' @param variance.component.list Optional precomputed list of variance components.
#'
#' @return The optimal lambda value.
#' @export
perform.cv = function(phenotype, genotype, input.kernel, MINQUE.type, KNN.type, lambda.list,
                      weight, fixed.effects, constrain, criteria) {

  sample.size = length(phenotype)
  n.lambda.list = length(lambda.list)

  # Split data for cross-validation
  batch.size = sample.size %/% n.lambda.list
  n.genotype.regions = length(genotype)
  randomized.indices = sample(1:sample.size)
  cv.fold.indices = list()

  for (fold.index in 0:(n.lambda.list - 2)) {
    fold.samples = randomized.indices[(fold.index * batch.size + 1):((fold.index + 1) * batch.size)]
    cv.fold.indices = c(cv.fold.indices, list(fold.samples))
  }
  last.fold.samples = randomized.indices[-c(1:((n.lambda.list - 1) * batch.size))]
  cv.fold.indices = c(cv.fold.indices, list(last.fold.samples))

  # Initialize evaluation metrics
  mse.values = numeric(n.lambda.list)
  correlation.values = numeric(n.lambda.list)
  loss.values = numeric(n.lambda.list)

  # Cross-validation loop
  for (lambda.index in 1:n.lambda.list) {
    train.indices = unlist(cv.fold.indices[-lambda.index])
    test.index = unlist(cv.fold.indices[lambda.index])
    phenotype.train = phenotype[train.indices]
    phenotype.test = phenotype[test.index]
    fixed.effects.train = fixed.effects[train.indices, , drop=FALSE]
    fixed.effects.test = fixed.effects[test.index, , drop=FALSE]
    genotype.train = list()
    genotype.test = list()
    genotype.test.loss = list()

    for (region.index in 1:n.genotype.regions) {
      genotype.train.region = genotype[[region.index]][train.indices, , drop=FALSE]
      genotype.test.region = genotype[[region.index]][test.index, , drop=FALSE]
      genotype.train = c(genotype.train, list(genotype.train.region))
      genotype.test = c(genotype.test, list(rbind(genotype.train.region, genotype.test.region)))
      genotype.test.loss = c(genotype.test.loss, list(genotype.test.region))
    }

    # Estimate model (using the fit.internal function to avoid infinite recursion)
    result = fit.internal(phenotype=phenotype.train, genotype=genotype.train, input.kernel=input.kernel,
                          MINQUE.type=MINQUE.type, KNN.type=KNN.type, lambda=lambda.list[lambda.index],
                          weight=weight, fixed.effects=fixed.effects.train, constrain=constrain)

    # Predict
    if (is.null(fixed.effects)) {
      phenotype.predicted = KNN.predict(phenotype=phenotype.train, genotype=genotype.test,
                                        input.kernel=input.kernel, theta=result$theta,
                                        KNN.type=KNN.type, beta=result$beta,
                                        fixed.effects.train=NULL, fixed.effects.test=NULL)
    } else {
      phenotype.predicted = KNN.predict(phenotype=phenotype.train, genotype=genotype.test,
                                        input.kernel=input.kernel, theta=result$theta,
                                        KNN.type=KNN.type, beta=result$beta,
                                        fixed.effects.train=fixed.effects.train,
                                        fixed.effects.test=fixed.effects.test)
    }

    # Calculate evaluation metrics
    mse.values[lambda.index] = mean((phenotype.test - phenotype.predicted)^2)
    correlation.values[lambda.index] = cor(phenotype.test, phenotype.predicted)
  }

  # Select optimal lambda based on criteria
  if (criteria == 'MSE') {
    optimal.index = which.min(mse.values)
  } else if (criteria == 'COR') {
    optimal.index = which.max(correlation.values)
  } else {
    stop("Unsupported criteria.")
  }
  lambda.optimal = lambda.list[optimal.index]
  if (is.na(lambda.optimal)) {
    lambda.optimal = 0
  }

  return(lambda.optimal)
}

############################################################
# Internal Fitting Function (to avoid recursion in CV)
############################################################

#' @title Internal Fitting Function
#'
#' @description
#' Internal function that performs the actual KNN fitting without CV.
#' This is used within the CV loop to avoid infinite recursion.
#'
fit.internal = function(phenotype, genotype, input.kernel, MINQUE.type, KNN.type, lambda,
                        weight, fixed.effects, constrain) {

  # If variance components are not provided, compute them
  variance.component.list = KNN2LMM(genotype=genotype, input.kernel=input.kernel, KNN.type=KNN.type)
  n.variance.component = length(variance.component.list)
  sample.size = length(phenotype)

  # Initialize V matrix based on MINQUE type
  if (MINQUE.type == 'MINQUE0') {
    variance.matrix = diag(1, sample.size, sample.size)
  } else if (MINQUE.type == 'MINQUE1') {
    if (is.null(weight)) {
      weight = rep(1 / n.variance.component, n.variance.component)
    }
    variance.matrix = Reduce(`+`, mapply(`*`, variance.component.list, weight, SIMPLIFY = FALSE))
  }

  # Ensure V matrix is positive definite
  npd.result = nearPD(variance.matrix)
  variance.matrix = as.matrix(npd.result$mat)
  precision.matrix = solve(variance.matrix)

  # Compute R*variance and R*phenotype
  rotated.variance.list = list()
  if (is.null(fixed.effects)) {
    beta = NULL
  } else {
    linear.model = lm(phenotype ~ fixed.effects - 1)
    beta = as.numeric(coef(linear.model))
    phenotype = residuals(linear.model)
  }
  for (component.index in seq.int(n.variance.component)) {
    rotated.variance.list[[component.index]] = precision.matrix %*% variance.component.list[[component.index]]
  }
  rotated.phenotype = precision.matrix %*% phenotype

  # Compute U matrix
  u.matrix = numeric(n.variance.component)
  for (component.index in seq.int(n.variance.component)) {
    u.matrix[component.index] = sum(rotated.phenotype * (variance.component.list[[component.index]] %*% rotated.phenotype))
  }

  # Compute C matrix
  covariance.matrix = matrix(0, n.variance.component, n.variance.component)
  for (row.index in seq.int(n.variance.component)) {
    for (col.index in seq.int(n.variance.component)) {
      covariance.matrix[row.index, col.index] = sum(rotated.variance.list[[row.index]] * rotated.variance.list[[col.index]])
    }
  }

  # Adjust C matrix with penalty
  penalty.matrix = diag(lambda, n.variance.component, n.variance.component)
  penalty.matrix[1, 1] = 0  # No penalty on the first variance component
  covariance.matrix = covariance.matrix + penalty.matrix

  # Ensure C matrix is positive definite
  npd.result = nearPD(covariance.matrix)
  covariance.matrix = as.matrix(npd.result$mat)

  # Estimate variance components
  covariance.matrix = solve(covariance.matrix)
  theta = covariance.matrix %*% u.matrix
  theta = as.numeric(theta)
  if (constrain) {
    theta[theta < 0] = 0
  }

  return(list(theta=theta, beta=beta))
}


############################################################
# KNN Numerical Estimation Function
############################################################

#' @title KNN Numerical Estimation
#'
#' @description
#' Estimates variance components in a KNN model using numerical optimization (e.g., glmnet).
#'
#' @param phenotype Numeric vector of response variable (phenotype) of length N.
#' @param genotype List of N x p genotype matrices.
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param MINQUE.type Type of MINQUE estimator ('MINQUE0' or 'MINQUE1') (default is 'MINQUE0').
#' @param KNN.type Type of KNN model ('KNN', 'LMM') (default is 'KNN').
#' @param lambda Regularization parameter.
#' @param weight Optional weight vector for variance components.
#' @param fixed.effects Optional fixed effects design matrix.
#' @param alpha Elastic net mixing parameter (default is 0).
#' @param lower.limits Lower limits for coefficients (default is 0).
#' @param variance.component.list Optional precomputed list of variance components.
#'
#' @return A list containing estimated variance components (`theta`), intercept (`intercept`), Q matrix (`precision.matrix`), and variance component list.


#' @title KNN Numerical Estimation
#'
#' @description
#' Estimates variance components in a KNN model via numerical optimization,
#' matching the closed-form KNN by whitening the data before fitting.
#' @export
KNN.numerical = function(phenotype,genotype,input.kernel = NULL,MINQUE.type = 'MINQUE0',
    KNN.type = 'KNN',lambda = NULL,weight = NULL,fixed.effects = NULL,
    alpha = 0,lower.limits = 0,variance.component.list = NULL) {
  # 1) Prepare variance components
  if (is.null(variance.component.list)) {
    variance.component.list = KNN2LMM(genotype = genotype,input.kernel = input.kernel,KNN.type = KNN.type)
  }
  n.variance.component = length(variance.component.list)
  sample.size = length(phenotype)

  # Initialize V matrix based on MINQUE type
  if (MINQUE.type == 'MINQUE0') {
    variance.matrix = diag(1, sample.size, sample.size)
  } else if (MINQUE.type == 'MINQUE1') {
    if (is.null(weight)) weight = rep(1 / n.variance.component, n.variance.component)
    variance.matrix = Reduce(`+`, mapply(`*`, variance.component.list, weight, SIMPLIFY = FALSE))
  }

  # Ensure positive definiteness
  npd.result = nearPD(variance.matrix)
  variance.matrix = as.matrix(npd.result$mat)

  # 3) Precision and whitening
  precision.matrix = solve(variance.matrix)
  precision.matrix.half = chol(precision.matrix)  # upper-triangular so that t(precision.matrix.half) %*% precision.matrix.half = precision.matrix

  # 4) Handle fixed effects if provided
  if (!is.null(fixed.effects)) {
    lm.fit    = lm(phenotype ~ fixed.effects - 1)
    beta      = as.numeric(coef(lm.fit))
    phenotype = residuals(lm.fit)
  } else {
    beta = NULL
  }

  # 5) Whiten phenotype and kernels
  rotated.phenotype = c(precision.matrix.half %*% phenotype %*% t(precision.matrix.half %*% phenotype))
  rotated.variance  = lapply(variance.component.list, function(kernel) {c(precision.matrix.half %*% kernel %*% t(precision.matrix.half))})
  design.matrix = do.call(cbind, rotated.variance)

  # 6) Penalty factors (no penalty on first component)
  penalty.factor = rep(1, n.variance.component)
  penalty.factor[1] = 0

  # 7) Cross-validation to choose lambda if needed
  cv.fit = cv.glmnet(design.matrix,rotated.phenotype,alpha = alpha,lower.limits = lower.limits,penalty.factor = penalty.factor)
  if (is.null(lambda)) {
    lambda.min = cv.fit$lambda.min
    sd = cv.fit$cvsd[which(cv.fit$lambda == lambda.min)]
    lambda.range = c(lambda.min - sd, lambda.min + sd)
    lambda.optimal = max(cv.fit$lambda[cv.fit$lambda >= lambda.range[1] & cv.fit$lambda <= lambda.range[2]])
  } else {
    lambda.optimal = lambda
  }

  # 8) Final model fit
  optimal.model = glmnet(design.matrix,rotated.phenotype,alpha = alpha,
    lower.limits = lower.limits,penalty.factor = penalty.factor,lambda = lambda.optimal)
  theta = as.numeric(coef(optimal.model))[-1]

  # 9) Return estimated components
  return(list(theta = theta,beta = beta,precision.matrix = precision.matrix, variance.component.list = variance.component.list))
}

############################################################
# KNN Prediction Function
############################################################

#' @title KNN Prediction
#'
#' @description
#' Predicts phenotypes for new subjects using the KNN model.
#'
#' @param phenotype Response variable (phenotype) of training subjects.
#' @param genotype List of N x p genotype matrices (combined training and testing subjects).
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param theta Estimated variance components.
#' @param KNN.type Type of KNN model ('KNN', 'LMM') (default is 'KNN').
#' @param beta Estimated fixed effects coefficients.
#' @param fixed.effects.train Fixed effects design matrix for training subjects.
#' @param fixed.effects.test Fixed effects design matrix for testing subjects.
#' @param variance.component.list Optional precomputed list of variance components.
#'
#' @return Predicted phenotypes for testing subjects.
#' @export
KNN.predict = function(phenotype, genotype, theta, input.kernel = NULL, KNN.type='KNN',
                       beta=NULL, fixed.effects.train=NULL, fixed.effects.test=NULL,
                       variance.component.list=NULL) {

  # Get variance components
  if (is.null(variance.component.list)) {
    variance.component.list = KNN2LMM(genotype=genotype, input.kernel=input.kernel, KNN.type=KNN.type)
  }
  n.variance.component = length(variance.component.list)

  total.sample.size = nrow(variance.component.list[[1]])
  training.sample.size = length(phenotype)
  test.sample.size = total.sample.size - training.sample.size

  # Initialize fixed effects if not provided
  if (is.null(fixed.effects.train) || is.null(fixed.effects.test)) {
    fixed.effects.train = matrix(0, training.sample.size, 1)
    fixed.effects.test = matrix(0, test.sample.size, 1)
    beta = 0
  }

  # Compute kernel matrix for random effects
  kernel.random.effects = Reduce(`+`, mapply(`*`, variance.component.list[-1], theta[-1], SIMPLIFY=FALSE))

  # Compute covariance matrices
  if (test.sample.size == 0) {
    covariance.test.train = kernel.random.effects[1:training.sample.size, 1:training.sample.size]
    covariance.train = kernel.random.effects[1:training.sample.size, 1:training.sample.size] +
      diag(theta[1], training.sample.size)
  } else {
    covariance.test.train = kernel.random.effects[(training.sample.size + 1):total.sample.size, 1:training.sample.size]
    covariance.train = kernel.random.effects[1:training.sample.size, 1:training.sample.size] +
      diag(theta[1], training.sample.size)
  }

  # Predict phenotype
  if (test.sample.size == 0) {
    phenotype.predicted = as.numeric(fixed.effects.train %*% beta +
                                       covariance.test.train %*% solve(covariance.train) %*%
                                       (phenotype - fixed.effects.train %*% beta))
  } else {
    phenotype.predicted = as.numeric(fixed.effects.test %*% beta +
                                       covariance.test.train %*% solve(covariance.train) %*%
                                       (phenotype - fixed.effects.train %*% beta))
  }

  return(phenotype.predicted)
}

############################################################
# KNN to LMM Conversion Function
############################################################

#' @title KNN to LMM Conversion
#'
#' @description
#' Converts kernel neural networks into linear mixed models by constructing variance component matrices.
#'
#' @param genotype List of N x p genotype matrices.
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param KNN.type Type of KNN model ('KNN', 'LMM').
#'
#' @return List of variance component matrices.
#' @export
KNN2LMM = function(genotype, input.kernel=NULL, KNN.type='KNN') {

  # Get the dimension of KNN
  sample.size = nrow(genotype[[1]])
  kernel.matrix.list = list()
  n.genotype.regions = length(genotype)

  # Generate kernel matrices
  if (is.null(input.kernel)) {
    input.kernel <- rep('product', n.genotype.regions)
  } else if (length(input.kernel) == 1) {
    input.kernel <- rep(input.kernel, n.genotype.regions)
  } else if (length(input.kernel) != n.genotype.regions) {
    stop(paste("input.kernel must have length 1 or", n.genotype.regions))
  }

  for (region.index in 1:n.genotype.regions) {
    kernel.matrix.temp <- findKernel(input.kernel[region.index], genotype[[region.index]])
    kernel.matrix.list <- c(kernel.matrix.list, list(kernel.matrix.temp))
  }
  n.kernel.matrices = length(kernel.matrix.list)

  # Initialize variance component list
  identity.matrix = list(diag(1, sample.size, sample.size))
  ones.matrix = list(matrix(1, sample.size, sample.size))
  variance.component.base = c(ones.matrix, kernel.matrix.list)
  variance.component.list = list()

  # Construct variance component list based on KNN.type
  if (KNN.type == 'KNN') {
    variance.component.list = c(identity.matrix, variance.component.base)
    # Add squared terms
    for (component.index in 2:(n.kernel.matrices + 1)){
      variance.component.squared = list(variance.component.base[[component.index]]^2)
      variance.component.list = c(variance.component.list, variance.component.squared)
    }
    # Add interaction terms
    if (n.kernel.matrices > 1){
      for (first.index in 2:(n.kernel.matrices + 1)){
        if (first.index < n.kernel.matrices + 1){
          for (second.index in (first.index + 1):(n.kernel.matrices + 1)){
            variance.component.interaction = list(variance.component.base[[first.index]] *
                                                    variance.component.base[[second.index]])
            variance.component.list = c(variance.component.list, variance.component.interaction)
          }
        }
      }
    }
  } else if (KNN.type == 'LMM') {
    variance.component.list = c(identity.matrix, kernel.matrix.list)
  }
  return(variance.component.list)
}

############################################################
# Generalized KNN (GKNN) Function with Cross-Validation
############################################################

#' @title Generalized Kernel Neural Network (GKNN) with Cross-Validation
#'
#' @description
#' Implements the GKNN model for non-Gaussian data types (e.g., binary, count data).
#' If lambda is not specified, performs cross-validation to select optimal lambda.
#'
#' @param phenotype Response variable (phenotype), a numeric vector of length N.
#' @param genotype List of N x p genotype matrices.
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param family The exponential family ('binomial', 'poisson', 'negbinomial', 'zinb') (default 'binomial').
#' @param MINQUE.type Type of MINQUE estimator ('MINQUE1' recommended for GKNN) (default 'MINQUE1').
#' @param KNN.type Type of KNN model ('KNN', 'LMM') (default 'KNN').
#' @param lambda Regularization parameter. If NULL, cross-validation will be performed.
#' @param lambda.list List of lambda values to try in cross-validation (used only if lambda is NULL).
#' @param weight Optional weight vector for variance components.
#' @param beta Initial estimates of fixed effects coefficients.
#' @param fixed.effects Fixed effects design matrix.
#' @param constrain Logical indicating whether to constrain variance components to be non-negative (default TRUE).
#' @param cv.criteria Criterion for selecting optimal lambda in CV ('DEVIANCE', 'AUC', 'MSE') (default is 'DEVIANCE').
#' @param thresh Convergence threshold for iterative methods (default 1e-3).
#' @param max.iter Maximum number of iterations (default 100).
#' @param alpha Elastic Net mixing parameter (used if numerical=TRUE) (default 0).
#' @param lower.limits Lower limits for variance components (used if numerical=TRUE) (default 0).
#' @param numerical Logical indicating whether to use numerical optimization methods (default FALSE).
#' @param theta.nb Overdispersion parameter for negative binomial (estimated if NULL).
#' @param pi.zinb Zero-inflation probability for ZINB (estimated if NULL).
#'
#' @return List containing estimated variance components, fixed effects, random effects, and other matrices.
#' @export
GKNN = function(phenotype, genotype, input.kernel = NULL, family='binomial', MINQUE.type='MINQUE1',
                KNN.type='KNN', lambda=NULL, lambda.list=NULL, weight=NULL, beta=NULL,
                fixed.effects=NULL, constrain=TRUE, cv.criteria='DEVIANCE', thresh=1e-3,
                max.iter=100, alpha=0, lower.limits=0, numerical=FALSE,
                theta.nb=NULL, pi.zinb=NULL) {

  sample.size = length(phenotype)

  # Check if fixed.effects is provided
  if (is.null(fixed.effects)) {
    stop("Design matrix fixed.effects must be provided.")
  }

  # If lambda is not specified, perform cross-validation
  if (is.null(lambda)) {
    # Initialize lambda list if not provided
    if (is.null(lambda.list)) {
      lambda.list = c(0.05, 5, 500, 50000, 5000000)
    }

    # Perform cross-validation for GKNN
    lambda.opt = perform.gknn.cv(phenotype=phenotype, genotype=genotype, input.kernel=input.kernel,
                                 family=family, MINQUE.type=MINQUE.type, KNN.type=KNN.type,
                                 lambda.list=lambda.list, weight=weight, beta=beta,
                                 fixed.effects=fixed.effects, constrain=constrain,
                                 criteria=cv.criteria, thresh=thresh, max.iter=max.iter,
                                 alpha=alpha, lower.limits=lower.limits, numerical=numerical,
                                 theta.nb=theta.nb, pi.zinb=pi.zinb)

    lambda = lambda.opt
  } else {
    lambda.opt = NULL
  }

  # Initialize beta if not provided
  if (is.null(beta)) {
    if (family %in% c('binomial', 'poisson')) {
      glm.result = glm(phenotype ~ fixed.effects - 1, family=family)
      beta.current = as.numeric(glm.result$coefficients)
    } else if (family == 'negbinomial') {
      # Use Poisson GLM for initial values
      glm.result = glm(phenotype ~ fixed.effects - 1, family='poisson')
      beta.current = as.numeric(glm.result$coefficients)
      # Initialize theta.nb if not provided
      if (is.null(theta.nb)) {
        theta.nb = 1  # Start with reasonable overdispersion
      }
    } else if (family == 'zinb') {
      # Use Poisson GLM for initial values
      glm.result = glm(phenotype[phenotype > 0] ~ fixed.effects[phenotype > 0, ] - 1, family='poisson')
      beta.current = as.numeric(glm.result$coefficients)
      # Initialize parameters if not provided
      if (is.null(theta.nb)) {
        theta.nb = 1  # Start with reasonable overdispersion
      }
      if (is.null(pi.zinb)) {
        pi.zinb = sum(phenotype == 0) / length(phenotype)  # Empirical zero proportion
      }
    } else {
      stop("Unsupported family.")
    }
  } else {
    beta.current = beta
  }

  random.phenotype = rnorm(sample.size)

  # Generate variance component list
  variance.component.list = KNN2LMM(genotype=genotype, input.kernel=input.kernel, KNN.type=KNN.type)
  variance.component.list.main = variance.component.list[-1]
  theta.current = rep(var(random.phenotype) / length(variance.component.list), length(variance.component.list))

  # Initialize working variables based on family
  if (family == 'binomial') {
    linear.predictor = fixed.effects %*% beta.current
    mean.response = exp(linear.predictor) / (1 + exp(linear.predictor))
    pseudo.phenotype = (phenotype - mean.response) / (mean.response * (1 - mean.response)) +
      log(mean.response / (1 - mean.response))
    working.weights = mean.response * (1 - mean.response)
  } else if (family == 'poisson') {
    linear.predictor = fixed.effects %*% beta.current
    mean.response = exp(linear.predictor)
    pseudo.phenotype = (phenotype - mean.response) / mean.response + log(mean.response)
    working.weights = mean.response
  } else if (family == 'negbinomial') {
    linear.predictor = fixed.effects %*% beta.current
    mean.response = exp(linear.predictor)
    variance.response = mean.response + mean.response^2 / theta.nb
    pseudo.phenotype = (phenotype - mean.response) / variance.response + log(mean.response)
    working.weights = mean.response / variance.response
  } else if (family == 'zinb') {
    linear.predictor = fixed.effects %*% beta.current
    mean.nb = exp(linear.predictor)
    # Expected value for ZINB: (1 - pi) * mu_nb
    mean.response = (1 - pi.zinb) * mean.nb
    # Variance for ZINB
    variance.response = (1 - pi.zinb) * (mean.nb + mean.nb^2 / theta.nb) +
      pi.zinb * (1 - pi.zinb) * mean.nb^2
    # Handle zero observations specially
    is.zero = (phenotype == 0)
    pseudo.phenotype = numeric(sample.size)
    working.weights = numeric(sample.size)
    # For non-zero observations
    pseudo.phenotype[!is.zero] = (phenotype[!is.zero] - mean.response[!is.zero]) /
      variance.response[!is.zero] + log(mean.nb[!is.zero])
    working.weights[!is.zero] = mean.response[!is.zero] / variance.response[!is.zero]
    # For zero observations
    pseudo.phenotype[is.zero] = log(mean.nb[is.zero])
    working.weights[is.zero] = mean.response[is.zero] / variance.response[is.zero]
  } else {
    stop("Unsupported family.")
  }

  # Iterative estimation loop
  for (iteration in 1:max.iter) {
    beta.previous = beta.current
    theta.previous = theta.current
    if (family == 'negbinomial') {
      theta.nb.previous = theta.nb
    }
    if (family == 'zinb') {
      theta.nb.previous = theta.nb
      pi.zinb.previous = pi.zinb
    }

    variance.component.list[[1]] = diag(as.numeric(1 / working.weights), sample.size, sample.size)

    # Estimate variance components
    if (numerical) {
      result = KNN.numerical(phenotype=random.phenotype, genotype=genotype, input.kernel=input.kernel,
                             MINQUE.type=MINQUE.type, KNN.type=KNN.type, lambda=lambda, weight=weight,
                             fixed.effects=NULL, alpha=alpha, lower.limits=lower.limits,
                             variance.component.list=variance.component.list)
    } else {
      result = KNN(phenotype=random.phenotype, genotype=genotype, input.kernel=input.kernel,
                   MINQUE.type=MINQUE.type, KNN.type=KNN.type, lambda=lambda, weight=weight,
                   fixed.effects=NULL, constrain=constrain, variance.component.list=variance.component.list)
    }
    theta.current = result$theta
    scale.parameter = theta.current[1]
    theta.random = theta.current[-1]

    # Update random effects covariance matrix
    random.effects.covariance = Reduce(`+`, mapply(`*`, variance.component.list.main, theta.random, SIMPLIFY=FALSE))
    npd.result = nearPD(scale.parameter * diag(as.numeric(1 / working.weights), sample.size, sample.size) + random.effects.covariance + diag(0.1, sample.size,sample.size))
    covariance.inverse = solve(as.matrix(npd.result$mat))

    # Update fixed effects beta
    beta.current = solve(t(fixed.effects) %*% covariance.inverse %*% fixed.effects) %*%
      t(fixed.effects) %*% covariance.inverse %*% pseudo.phenotype

    # Update random effects
    random.phenotype = random.effects.covariance %*% covariance.inverse %*%
      (pseudo.phenotype - fixed.effects %*% beta.current)

    # Update parameters specific to NB and ZINB
    if (family == 'negbinomial') {
      # Update theta.nb using method of moments
      linear.predictor = fixed.effects %*% beta.current + random.phenotype
      mean.response = exp(linear.predictor)
      residuals.sq = (phenotype - mean.response)^2
      theta.nb = mean(mean.response^2 / pmax(residuals.sq - mean.response, 0.1))
      theta.nb = pmax(theta.nb, 0.1)  # Ensure positive
    } else if (family == 'zinb') {
      # Update theta.nb and pi.zinb
      linear.predictor = fixed.effects %*% beta.current + random.phenotype
      mean.nb = exp(linear.predictor)

      # EM-like update for pi.zinb
      is.zero = (phenotype == 0)
      prob.zero.from.zinb = pi.zinb / (pi.zinb + (1 - pi.zinb) * dnbinom(0, size=theta.nb, mu=mean.nb))
      pi.zinb = mean(is.zero * prob.zero.from.zinb + (1 - is.zero) * 0)
      pi.zinb = pmax(pmin(pi.zinb, 0.99), 0.01)  # Keep in (0.01, 0.99)

      # Update theta.nb for non-zero observations
      if (sum(!is.zero) > 0) {
        mean.response.nz = mean.nb[!is.zero]
        residuals.sq.nz = (phenotype[!is.zero] - mean.response.nz)^2
        theta.nb = mean(mean.response.nz^2 / pmax(residuals.sq.nz - mean.response.nz, 0.1))
        theta.nb = pmax(theta.nb, 0.1)  # Ensure positive
      }
    }

    # Update working variables based on family
    if (family == 'binomial') {
      linear.predictor = fixed.effects %*% beta.current + random.phenotype
      mean.response = exp(linear.predictor) / (1 + exp(linear.predictor))
      pseudo.phenotype = (phenotype - mean.response) / (mean.response * (1 - mean.response)) +
        log(mean.response / (1 - mean.response))
      working.weights = mean.response * (1 - mean.response)
    } else if (family == 'poisson') {
      linear.predictor = fixed.effects %*% beta.current + random.phenotype
      mean.response = exp(linear.predictor)
      pseudo.phenotype = (phenotype - mean.response) / mean.response + log(mean.response)
      working.weights = mean.response
    } else if (family == 'negbinomial') {
      linear.predictor = fixed.effects %*% beta.current + random.phenotype
      mean.response = exp(linear.predictor)
      variance.response = mean.response + mean.response^2 / theta.nb
      pseudo.phenotype = (phenotype - mean.response) / variance.response + log(mean.response)
      working.weights = mean.response / variance.response
    } else if (family == 'zinb') {
      linear.predictor = fixed.effects %*% beta.current + random.phenotype
      mean.nb = exp(linear.predictor)
      mean.response = (1 - pi.zinb) * mean.nb
      variance.response = (1 - pi.zinb) * (mean.nb + mean.nb^2 / theta.nb) +
        pi.zinb * (1 - pi.zinb) * mean.nb^2
      is.zero = (phenotype == 0)
      pseudo.phenotype[!is.zero] = (phenotype[!is.zero] - mean.response[!is.zero]) /
        variance.response[!is.zero] + log(mean.nb[!is.zero])
      working.weights[!is.zero] = mean.response[!is.zero] / variance.response[!is.zero]
      pseudo.phenotype[is.zero] = log(mean.nb[is.zero])
      working.weights[is.zero] = mean.response[is.zero] / variance.response[is.zero]
    }

    # Check convergence every 10 iterations
    if (iteration %% 10 == 0) {
      convergence.difference = sum(abs(beta.current - beta.previous)) + sum(abs(theta.current - theta.previous))
      if (family == 'negbinomial') {
        convergence.difference = convergence.difference + abs(theta.nb - theta.nb.previous)
      }
      if (family == 'zinb') {
        convergence.difference = convergence.difference + abs(theta.nb - theta.nb.previous) + abs(pi.zinb - pi.zinb.previous)
      }
      if (convergence.difference < thresh) break
    }
  }

  theta = as.numeric(theta.current)
  beta = as.numeric(beta.current)
  random.phenotype = as.numeric(random.phenotype)
  if (constrain) theta[theta < 0] = 0

  result = list(theta=theta, beta=beta, random.phenotype=random.phenotype, covariance.inverse=covariance.inverse)

  if (family == 'negbinomial') {
    result$theta.nb = theta.nb
  }
  if (family == 'zinb') {
    result$theta.nb = theta.nb
    result$pi.zinb = pi.zinb
  }

  return(result)
}

############################################################
# GKNN Prediction Function
############################################################

#' @title GKNN Prediction
#'
#' @description
#' Predicts phenotypes for new subjects using the GKNN model.
#'
#' @param phenotype Response variable (phenotype) of training subjects.
#' @param genotype List of N x p genotype matrices (combined training and testing subjects).
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param theta Estimated variance components.
#' @param random.phenotype Estimated random effects.
#' @param covariance.inverse Estimated inverse covariance matrix of random effects.
#' @param family The exponential family ('binomial', 'poisson', 'negbinomial', 'zinb') (default 'binomial').
#' @param KNN.type Type of KNN model ('KNN', 'LMM') (default 'KNN').
#' @param beta Estimated fixed effects coefficients.
#' @param fixed.effects.train Fixed effects design matrix for training subjects.
#' @param fixed.effects.test Fixed effects design matrix for testing subjects.
#' @param theta.nb Overdispersion parameter for negative binomial (used if family is 'negbinomial' or 'zinb').
#' @param pi.zinb Zero-inflation probability for ZINB (used if family is 'zinb').
#'
#' @return Predicted phenotypes for testing subjects.
#' @export
GKNN.predict = function(phenotype, genotype, theta, random.phenotype, covariance.inverse,
                        input.kernel = NULL, family='binomial', KNN.type='KNN', beta=NULL,
                        fixed.effects.train=NULL, fixed.effects.test=NULL,
                        theta.nb=NULL, pi.zinb=NULL) {

  total.sample.size = nrow(genotype[[1]])
  training.sample.size = length(phenotype)
  test.sample.size = total.sample.size - training.sample.size

  # Initialize fixed effects if not provided
  if (is.null(fixed.effects.train) || is.null(fixed.effects.test)) {
    fixed.effects.train = matrix(0, training.sample.size, 1)
    fixed.effects.test = matrix(0, test.sample.size, 1)
    beta = 0
  }

  # Generate variance component list
  variance.component.list = KNN2LMM(genotype=genotype, input.kernel=input.kernel, KNN.type=KNN.type)
  n.variance.component = length(variance.component.list)

  # Compute kernel matrix for random effects
  kernel.random.effects = Reduce(`+`, mapply(`*`, variance.component.list[-1], theta[-1], SIMPLIFY=FALSE))

  # Compute covariance matrix between test and train
  if (test.sample.size == 0) {
    covariance.test.train = kernel.random.effects[1:training.sample.size, 1:training.sample.size]
  } else {
    covariance.test.train = kernel.random.effects[(training.sample.size + 1):total.sample.size, 1:training.sample.size]
  }

  # Predict random effects
  random.phenotype.predicted = covariance.test.train %*% covariance.inverse %*% random.phenotype

  # Compute predictions based on family
  if (family == 'binomial') {
    if (test.sample.size == 0) {
      linear.predictor = fixed.effects.train %*% beta + random.phenotype.predicted
      phenotype.predicted = as.numeric(exp(linear.predictor) / (1 + exp(linear.predictor)))
    } else {
      linear.predictor = fixed.effects.test %*% beta + random.phenotype.predicted
      phenotype.predicted = as.numeric(exp(linear.predictor) / (1 + exp(linear.predictor)))
    }
  } else if (family == 'poisson') {
    if (test.sample.size == 0) {
      linear.predictor = fixed.effects.train %*% beta + random.phenotype.predicted
      phenotype.predicted = as.numeric(exp(linear.predictor))
    } else {
      linear.predictor = fixed.effects.test %*% beta + random.phenotype.predicted
      phenotype.predicted = as.numeric(exp(linear.predictor))
    }
  } else if (family == 'negbinomial') {
    # For negative binomial, return the mean of the NB distribution
    if (test.sample.size == 0) {
      linear.predictor = fixed.effects.train %*% beta + random.phenotype.predicted
      phenotype.predicted = as.numeric(exp(linear.predictor))
    } else {
      linear.predictor = fixed.effects.test %*% beta + random.phenotype.predicted
      phenotype.predicted = as.numeric(exp(linear.predictor))
    }
  } else if (family == 'zinb') {
    # For ZINB, return the expected value: (1 - pi) * mu_nb
    if (test.sample.size == 0) {
      linear.predictor = fixed.effects.train %*% beta + random.phenotype.predicted
      mean.nb = as.numeric(exp(linear.predictor))
      phenotype.predicted = (1 - pi.zinb) * mean.nb
    } else {
      linear.predictor = fixed.effects.test %*% beta + random.phenotype.predicted
      mean.nb = as.numeric(exp(linear.predictor))
      phenotype.predicted = (1 - pi.zinb) * mean.nb
    }
  } else {
    stop("Unsupported family.")
  }
  return(phenotype.predicted)
}

############################################################
# GKNN Cross-Validation Function
############################################################

#' @title GKNN Cross-Validation Function
#'
#' @description
#' Performs cross-validation for lambda selection in GKNN models.
#' @export
perform.gknn.cv = function(phenotype, genotype, input.kernel, family, MINQUE.type, KNN.type,
                           lambda.list, weight, beta, fixed.effects, constrain, criteria,
                           thresh, max.iter, alpha, lower.limits, numerical,
                           theta.nb=NULL, pi.zinb=NULL) {

  sample.size = length(phenotype)
  n.lambda.list = length(lambda.list)

  # Split data for cross-validation
  batch.size = sample.size %/% n.lambda.list
  n.genotype.regions = length(genotype)
  randomized.indices = sample(1:sample.size)
  cv.fold.indices = list()

  for (fold.index in 0:(n.lambda.list - 2)) {
    fold.samples = randomized.indices[(fold.index * batch.size + 1):((fold.index + 1) * batch.size)]
    cv.fold.indices = c(cv.fold.indices, list(fold.samples))
  }
  last.fold.samples = randomized.indices[-c(1:((n.lambda.list - 1) * batch.size))]
  cv.fold.indices = c(cv.fold.indices, list(last.fold.samples))

  # Initialize evaluation metrics
  deviance.values = numeric(n.lambda.list)
  auc.values = numeric(n.lambda.list)
  mse.values = numeric(n.lambda.list)

  # Cross-validation loop
  for (lambda.index in 1:n.lambda.list) {
    train.indices = unlist(cv.fold.indices[-lambda.index])
    test.index = unlist(cv.fold.indices[lambda.index])
    phenotype.train = phenotype[train.indices]
    phenotype.test = phenotype[test.index]
    fixed.effects.train = fixed.effects[train.indices, , drop=FALSE]
    fixed.effects.test = fixed.effects[test.index, , drop=FALSE]
    genotype.train = list()
    genotype.test = list()

    for (region.index in 1:n.genotype.regions) {
      genotype.train.region = genotype[[region.index]][train.indices, , drop=FALSE]
      genotype.test.region = genotype[[region.index]][test.index, , drop=FALSE]
      genotype.train = c(genotype.train, list(genotype.train.region))
      genotype.test = c(genotype.test, list(rbind(genotype.train.region, genotype.test.region)))
    }

    # Fit GKNN model with internal function to avoid recursion
    result = fit.gknn.internal(phenotype=phenotype.train, genotype=genotype.train,
                               input.kernel=input.kernel, family=family, MINQUE.type=MINQUE.type,
                               KNN.type=KNN.type, lambda=lambda.list[lambda.index], weight=weight,
                               beta=beta, fixed.effects=fixed.effects.train, constrain=constrain,
                               thresh=thresh, max.iter=max.iter, alpha=alpha,
                               lower.limits=lower.limits, numerical=numerical,
                               theta.nb=theta.nb, pi.zinb=pi.zinb)

    # Predict
    phenotype.predicted = GKNN.predict(phenotype=phenotype.train, genotype=genotype.test,
                                       theta=result$theta, random.phenotype=result$random.phenotype,
                                       covariance.inverse=result$covariance.inverse,
                                       input.kernel=input.kernel, family=family,
                                       KNN.type=KNN.type, beta=result$beta,
                                       fixed.effects.train=fixed.effects.train,
                                       fixed.effects.test=fixed.effects.test,
                                       theta.nb=result$theta.nb, pi.zinb=result$pi.zinb)

    # Calculate evaluation metrics
    if (family == 'binomial') {
      # Deviance
      deviance.values[lambda.index] = -2 * sum(phenotype.test * log(phenotype.predicted + 1e-8) +
                                                 (1 - phenotype.test) * log(1 - phenotype.predicted + 1e-8))
      # AUC (simple approximation)
      if (criteria == 'AUC' && length(unique(phenotype.test)) == 2) {
        positive.predictions = phenotype.predicted[phenotype.test == 1]
        negative.predictions = phenotype.predicted[phenotype.test == 0]
        auc.values[lambda.index] = mean(outer(positive.predictions, negative.predictions, ">"))
      }
    } else if (family == 'poisson') {
      # Deviance
      deviance.values[lambda.index] = 2 * sum(phenotype.test * log((phenotype.test + 1e-8) /
                                                                     (phenotype.predicted + 1e-8)) -
                                                (phenotype.test - phenotype.predicted))
    } else if (family == 'negbinomial') {
      # Negative binomial deviance
      deviance.values[lambda.index] = 2 * sum(phenotype.test * log((phenotype.test + 1e-8) /
                                                                     (phenotype.predicted + 1e-8)) -
                                                (phenotype.test + result$theta.nb) *
                                                log((phenotype.test + result$theta.nb) /
                                                      (phenotype.predicted + result$theta.nb)))
    } else if (family == 'zinb') {
      # ZINB deviance (simplified)
      deviance.values[lambda.index] = -2 * sum(ifelse(phenotype.test == 0,
                                                      log(result$pi.zinb + (1 - result$pi.zinb) *
                                                            dnbinom(0, size=result$theta.nb, mu=phenotype.predicted)),
                                                      log((1 - result$pi.zinb) *
                                                            dnbinom(phenotype.test, size=result$theta.nb, mu=phenotype.predicted))))
    }
    mse.values[lambda.index] = mean((phenotype.test - phenotype.predicted)^2)
  }

  # Select optimal lambda based on criteria
  if (criteria == 'DEVIANCE') {
    optimal.index = which.min(deviance.values)
  } else if (criteria == 'AUC') {
    optimal.index = which.max(auc.values)
  } else if (criteria == 'MSE') {
    optimal.index = which.min(mse.values)
  } else {
    stop("Unsupported criteria for GKNN.")
  }

  lambda.optimal = lambda.list[optimal.index]
  if (is.na(lambda.optimal)) {
    lambda.optimal = 0
  }

  return(lambda.optimal)
}

############################################################
# Internal GKNN Fitting Function (to avoid recursion in CV)
############################################################

fit.gknn.internal = function(phenotype, genotype, input.kernel, family, MINQUE.type, KNN.type,
                             lambda, weight, beta, fixed.effects, constrain, thresh, max.iter,
                             alpha, lower.limits, numerical, theta.nb=NULL, pi.zinb=NULL) {

  # This is essentially the same as GKNN but without CV logic
  sample.size = length(phenotype)

  # Initialize beta if not provided
  if (is.null(beta)) {
    if (family %in% c('binomial', 'poisson')) {
      glm.result = glm(phenotype ~ fixed.effects - 1, family=family)
      beta.current = as.numeric(glm.result$coefficients)
    } else if (family == 'negbinomial') {
      glm.result = glm(phenotype ~ fixed.effects - 1, family='poisson')
      beta.current = as.numeric(glm.result$coefficients)
      if (is.null(theta.nb)) {
        theta.nb = 1
      }
    } else if (family == 'zinb') {
      glm.result = glm(phenotype[phenotype > 0] ~ fixed.effects[phenotype > 0, ] - 1, family='poisson')
      beta.current = as.numeric(glm.result$coefficients)
      if (is.null(theta.nb)) {
        theta.nb = 1
      }
      if (is.null(pi.zinb)) {
        pi.zinb = sum(phenotype == 0) / length(phenotype)
      }
    }
  } else {
    beta.current = beta
  }

  random.phenotype = rnorm(sample.size)

  # Generate variance component list
  variance.component.list = KNN2LMM(genotype=genotype, input.kernel=input.kernel, KNN.type=KNN.type)
  variance.component.list.main = variance.component.list[-1]
  theta.current = rep(var(random.phenotype) / length(variance.component.list), length(variance.component.list))

  # Initialize working variables based on family
  if (family == 'binomial') {
    linear.predictor = fixed.effects %*% beta.current
    mean.response = exp(linear.predictor) / (1 + exp(linear.predictor))
    pseudo.phenotype = (phenotype - mean.response) / (mean.response * (1 - mean.response)) +
      log(mean.response / (1 - mean.response))
    working.weights = mean.response * (1 - mean.response)
  } else if (family == 'poisson') {
    linear.predictor = fixed.effects %*% beta.current
    mean.response = exp(linear.predictor)
    pseudo.phenotype = (phenotype - mean.response) / mean.response + log(mean.response)
    working.weights = mean.response
  } else if (family == 'negbinomial') {
    linear.predictor = fixed.effects %*% beta.current
    mean.response = exp(linear.predictor)
    variance.response = mean.response + mean.response^2 / theta.nb
    pseudo.phenotype = (phenotype - mean.response) / variance.response + log(mean.response)
    working.weights = mean.response / variance.response
  } else if (family == 'zinb') {
    linear.predictor = fixed.effects %*% beta.current
    mean.nb = exp(linear.predictor)
    mean.response = (1 - pi.zinb) * mean.nb
    variance.response = (1 - pi.zinb) * (mean.nb + mean.nb^2 / theta.nb) +
      pi.zinb * (1 - pi.zinb) * mean.nb^2
    is.zero = (phenotype == 0)
    pseudo.phenotype = numeric(sample.size)
    working.weights = numeric(sample.size)
    pseudo.phenotype[!is.zero] = (phenotype[!is.zero] - mean.response[!is.zero]) /
      variance.response[!is.zero] + log(mean.nb[!is.zero])
    working.weights[!is.zero] = mean.response[!is.zero] / variance.response[!is.zero]
    pseudo.phenotype[is.zero] = log(mean.nb[is.zero])
    working.weights[is.zero] = mean.response[is.zero] / variance.response[is.zero]
  }

  # Iterative estimation loop
  for (iteration in 1:max.iter) {
    beta.previous = beta.current
    theta.previous = theta.current
    if (family == 'negbinomial') {
      theta.nb.previous = theta.nb
    }
    if (family == 'zinb') {
      theta.nb.previous = theta.nb
      pi.zinb.previous = pi.zinb
    }

    variance.component.list[[1]] = diag(as.numeric(1 / working.weights),sample.size,sample.size)

    # Estimate variance components
    if (numerical) {
      result = KNN.numerical(phenotype=random.phenotype, genotype=genotype, input.kernel=input.kernel,
                             MINQUE.type=MINQUE.type, KNN.type=KNN.type, lambda=lambda, weight=weight,
                             fixed.effects=NULL, alpha=alpha, lower.limits=lower.limits,
                             variance.component.list=variance.component.list)
    } else {
      result = KNN(phenotype=random.phenotype, genotype=genotype, input.kernel=input.kernel,
                   MINQUE.type=MINQUE.type, KNN.type=KNN.type, lambda=lambda, weight=weight,
                   fixed.effects=NULL, constrain=constrain,
                   variance.component.list=variance.component.list)
    }
    theta.current = result$theta
    scale.parameter = theta.current[1]
    theta.random = theta.current[-1]

    # Update random effects covariance matrix
    random.effects.covariance = Reduce(`+`, mapply(`*`, variance.component.list.main, theta.random, SIMPLIFY=FALSE))
    npd.result = nearPD(scale.parameter * diag(as.numeric(1 / working.weights), sample.size, sample.size) + random.effects.covariance)
    covariance.inverse = solve(as.matrix(npd.result$mat))

    # Update fixed effects beta
    beta.current = solve(t(fixed.effects) %*% covariance.inverse %*% fixed.effects) %*%
      t(fixed.effects) %*% covariance.inverse %*% pseudo.phenotype

    # Update random effects
    random.phenotype = random.effects.covariance %*% covariance.inverse %*%
      (pseudo.phenotype - fixed.effects %*% beta.current)

    # Update parameters specific to NB and ZINB
    if (family == 'negbinomial') {
      # Update theta.nb using method of moments
      linear.predictor = fixed.effects %*% beta.current + random.phenotype
      mean.response = exp(linear.predictor)
      residuals.sq = (phenotype - mean.response)^2
      theta.nb = mean(mean.response^2 / pmax(residuals.sq - mean.response, 0.1))
      theta.nb = pmax(theta.nb, 0.1)  # Ensure positive
    } else if (family == 'zinb') {
      # Update theta.nb and pi.zinb
      linear.predictor = fixed.effects %*% beta.current + random.phenotype
      mean.nb = exp(linear.predictor)

      # EM-like update for pi.zinb
      is.zero = (phenotype == 0)
      prob.zero.from.zinb = pi.zinb / (pi.zinb + (1 - pi.zinb) * dnbinom(0, size=theta.nb, mu=mean.nb))
      pi.zinb = mean(is.zero * prob.zero.from.zinb + (1 - is.zero) * 0)
      pi.zinb = pmax(pmin(pi.zinb, 0.99), 0.01)  # Keep in (0.01, 0.99)

      # Update theta.nb for non-zero observations
      if (sum(!is.zero) > 0) {
        mean.response.nz = mean.nb[!is.zero]
        residuals.sq.nz = (phenotype[!is.zero] - mean.response.nz)^2
        theta.nb = mean(mean.response.nz^2 / pmax(residuals.sq.nz - mean.response.nz, 0.1))
        theta.nb = pmax(theta.nb, 0.1)  # Ensure positive
      }
    }

    # Update working variables based on family
    if (family == 'binomial') {
      linear.predictor = fixed.effects %*% beta.current + random.phenotype
      mean.response = exp(linear.predictor) / (1 + exp(linear.predictor))
      pseudo.phenotype = (phenotype - mean.response) / (mean.response * (1 - mean.response)) +
        log(mean.response / (1 - mean.response))
      working.weights = mean.response * (1 - mean.response)
    } else if (family == 'poisson') {
      linear.predictor = fixed.effects %*% beta.current + random.phenotype
      mean.response = exp(linear.predictor)
      pseudo.phenotype = (phenotype - mean.response) / mean.response + log(mean.response)
      working.weights = mean.response
    } else if (family == 'negbinomial') {
      linear.predictor = fixed.effects %*% beta.current + random.phenotype
      mean.response = exp(linear.predictor)
      variance.response = mean.response + mean.response^2 / theta.nb
      pseudo.phenotype = (phenotype - mean.response) / variance.response + log(mean.response)
      working.weights = mean.response / variance.response
    } else if (family == 'zinb') {
      linear.predictor = fixed.effects %*% beta.current + random.phenotype
      mean.nb = exp(linear.predictor)
      mean.response = (1 - pi.zinb) * mean.nb
      variance.response = (1 - pi.zinb) * (mean.nb + mean.nb^2 / theta.nb) +
        pi.zinb * (1 - pi.zinb) * mean.nb^2
      is.zero = (phenotype == 0)
      pseudo.phenotype[!is.zero] = (phenotype[!is.zero] - mean.response[!is.zero]) /
        variance.response[!is.zero] + log(mean.nb[!is.zero])
      working.weights[!is.zero] = mean.response[!is.zero] / variance.response[!is.zero]
      pseudo.phenotype[is.zero] = log(mean.nb[is.zero])
      working.weights[is.zero] = mean.response[is.zero] / variance.response[is.zero]
    }

    # Check convergence every 10 iterations
    if (iteration %% 10 == 0) {
      convergence.difference = sum(abs(beta.current - beta.previous)) + sum(abs(theta.current - theta.previous))
      if (family == 'negbinomial') {
        convergence.difference = convergence.difference + abs(theta.nb - theta.nb.previous)
      }
      if (family == 'zinb') {
        convergence.difference = convergence.difference + abs(theta.nb - theta.nb.previous) + abs(pi.zinb - pi.zinb.previous)
      }
      if (convergence.difference < thresh) break
    }
  }

  theta = as.numeric(theta.current)
  beta = as.numeric(beta.current)
  random.phenotype = as.numeric(random.phenotype)
  if (constrain) theta[theta < 0] = 0

  result = list(theta=theta, beta=beta, random.phenotype=random.phenotype, covariance.inverse=covariance.inverse)

  if (family == 'negbinomial') {
    result$theta.nb = theta.nb
  }
  if (family == 'zinb') {
    result$theta.nb = theta.nb
    result$pi.zinb = pi.zinb
  }

  return(result)
}

############################################################
# KNN Cross Estimation Function
############################################################

#' @title KNN Cross Estimation
#'
#' @description
#' Estimates the cross-covariance between two phenotypes using KNN.
#'
#' @param phenotype.first Numeric vector of the first phenotype, length N.
#' @param phenotype.second Numeric vector of the second phenotype, length N.
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param variance.component.list List of variance components.
#' @param MINQUE.type Type of MINQUE estimator ('MINQUE0' or 'MINQUE1') (default 'MINQUE0').
#' @param lambda Regularization parameter (default 0).
#' @param weight Optional weight vector for variance components.
#' @param constrain Logical indicating whether to constrain variance components to be non-negative (default FALSE).
#'
#' @return A list containing estimated variance components (`theta`), C matrix (`covariance.matrix`), Q matrix (`precision.matrix`), and variance component list.
#' @export
KNN.cross = function(phenotype.first, phenotype.second, input.kernel = NULL, variance.component.list,
                     MINQUE.type='MINQUE0', lambda=0, weight=NULL, constrain=FALSE) {

  n.variance.component = length(variance.component.list)
  sample.size = length(phenotype.first)

  # Initialize V matrix based on MINQUE type
  if (MINQUE.type == 'MINQUE0') {
    variance.matrix = diag(1, sample.size, sample.size)
  } else if (MINQUE.type == 'MINQUE1') {
    if (is.null(weight)) {
      weight = rep(1 / n.variance.component, n.variance.component)
    }
    variance.matrix = Reduce(`+`, mapply(`*`, variance.component.list, weight, SIMPLIFY=FALSE))
  }

  # Ensure V matrix is positive definite
  npd.result = nearPD(variance.matrix)
  variance.matrix = as.matrix(npd.result$mat)
  precision.matrix = solve(variance.matrix)

  # Compute R*variance and R*phenotype
  rotated.variance.list = list()
  for (component.index in seq_len(n.variance.component)) {
    rotated.variance.list[[component.index]] = precision.matrix %*% variance.component.list[[component.index]]
  }
  rotated.phenotype.first = precision.matrix %*% phenotype.first
  rotated.phenotype.second = precision.matrix %*% phenotype.second

  # Compute U matrix (cross-product)
  u.matrix = numeric(n.variance.component)
  for (component.index in seq_len(n.variance.component)) {
    u.matrix[component.index] = sum(rotated.phenotype.first *
                                      (variance.component.list[[component.index]] %*% rotated.phenotype.second))
  }

  # Compute C matrix
  covariance.matrix = matrix(0, n.variance.component, n.variance.component)
  for (row.index in seq_len(n.variance.component)) {
    for (col.index in seq_len(n.variance.component)) {
      covariance.matrix[row.index, col.index] = sum(rotated.variance.list[[row.index]] *
                                                      rotated.variance.list[[col.index]])
    }
  }

  # Apply penalty to C matrix
  penalty.matrix = diag(lambda, n.variance.component)
  penalty.matrix[1, 1] = 0  # No penalty on the first variance component
  covariance.matrix = covariance.matrix + penalty.matrix
  npd.result = nearPD(covariance.matrix)
  covariance.matrix = as.matrix(npd.result$mat)

  # Estimate variance components
  covariance.matrix = solve(covariance.matrix)
  theta = covariance.matrix %*% u.matrix
  theta = as.numeric(theta)
  if (constrain) {
    theta[theta < 0] = 0
  }

  return(list(theta=theta, covariance.matrix=covariance.matrix, precision.matrix=precision.matrix,
              variance.component.list=variance.component.list))
}

############################################################
# Multivariate KNN (MVKNN) Prediction Function
############################################################

#' @title MVKNN Prediction
#'
#' @description
#' Predicts multivariate phenotypes for new subjects using the MVKNN model.
#'
#' @param phenotype Matrix of training phenotypes (N.train x d), where N.train is the number of training samples and d is the number of traits.
#' @param genotype List of N x p genotype matrices (combined training and testing subjects).
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param theta Estimated variance components.
#' @param phenotype.covariance Estimated covariance matrix of phenotypes (d x d).
#' @param KNN.type Type of KNN model ('KNN', 'LMM') (default 'KNN').
#' @param beta Estimated fixed effects coefficients (optional).
#' @param fixed.effects.train Fixed effects design matrix for training subjects (optional).
#' @param fixed.effects.test Fixed effects design matrix for testing subjects (optional).
#'
#' @return Predicted phenotypes for testing subjects as a matrix (N.test x d).
#' @export
MVKNN.predict = function(phenotype, genotype, theta, phenotype.covariance, input.kernel = NULL,
                         KNN.type='KNN', beta=NULL, fixed.effects.train=NULL, fixed.effects.test=NULL) {

  total.sample.size = nrow(genotype[[1]])
  training.sample.size = nrow(phenotype)
  test.sample.size = total.sample.size - training.sample.size
  n.traits = ncol(phenotype)
  total.dimension = total.sample.size * n.traits
  training.dimension = training.sample.size * n.traits
  test.dimension = test.sample.size * n.traits
  identity.sample = diag(total.sample.size)
  identity.trait = diag(n.traits)
  phenotype.vectorized = as.vector(t(phenotype))  # Convert phenotype to a vector

  # Initialize fixed effects if not provided
  if (is.null(fixed.effects.train) || is.null(fixed.effects.test)) {
    fixed.effects.train = matrix(0, nrow=training.sample.size, ncol=1)
    fixed.effects.test = matrix(0, nrow=test.sample.size, ncol=1)
    beta = matrix(0, nrow=1, ncol=n.traits)
  }

  # Get variance components
  variance.component.list = KNN2LMM(genotype=genotype, input.kernel=input.kernel, KNN.type=KNN.type)
  n.variance.component = length(variance.component.list)

  # Adjust variance components for multivariate case
  variance.component.list[[1]] = kronecker(identity.sample, phenotype.covariance)
  for (component.index in 2:n.variance.component) {
    variance.component.list[[component.index]] = kronecker(variance.component.list[[component.index]],
                                                           phenotype.covariance)
  }

  # Compute kernel matrix for random effects
  kernel.random.effects = Reduce(`+`, mapply(`*`, variance.component.list[-1], theta[-1], SIMPLIFY=FALSE))

  # Compute covariance matrices
  if (test.sample.size == 0) {
    covariance.test.train = kernel.random.effects[1:training.dimension, 1:training.dimension]
    covariance.train = kernel.random.effects[1:training.dimension, 1:training.dimension] +
      kronecker(diag(training.sample.size), phenotype.covariance)
  } else {
    covariance.test.train = kernel.random.effects[(training.dimension + 1):total.dimension, 1:training.dimension]
    covariance.train = kernel.random.effects[1:training.dimension, 1:training.dimension] +
      kronecker(diag(training.sample.size), phenotype.covariance)
  }

  # Predict phenotype
  if (test.sample.size == 0) {
    fixed.effect.vector = as.vector(t(beta) %*% t(fixed.effects.train))
    phenotype.predicted.vector = fixed.effect.vector +
      covariance.test.train %*% solve(covariance.train) %*%
      (phenotype.vectorized - fixed.effect.vector)
  } else {
    fixed.effect.train.vector = as.vector(t(beta) %*% t(fixed.effects.train))
    fixed.effect.test.vector = as.vector(t(beta) %*% t(fixed.effects.test))
    phenotype.predicted.vector = fixed.effect.test.vector +
      covariance.test.train %*% solve(covariance.train) %*%
      (phenotype.vectorized - fixed.effect.train.vector)
  }

  # Convert predicted phenotype vector back to matrix form
  phenotype.predicted = t(matrix(phenotype.predicted.vector, nrow=n.traits, ncol=test.sample.size))

  # Return predicted phenotypes for testing subjects
  return(phenotype.predicted)
}

############################################################
# Multivariate KNN (MVKNN) Estimation Function with CV
############################################################

#' @title Multivariate KNN (MVKNN) Estimation with Cross-Validation
#'
#' @description
#' Estimates variance components and phenotype covariance in a multivariate KNN model.
#' If lambda is not specified, performs cross-validation to select optimal lambda.
#'
#' @param phenotype Matrix of phenotypes (N x d), where N is the number of samples and d is the number of traits.
#' @param genotype List of N x p genotype matrices.
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param MINQUE.type Type of MINQUE estimator ('MINQUE0' or 'MINQUE1') (default 'MINQUE0').
#' @param KNN.type Type of KNN model ('KNN', 'LMM') (default 'KNN').
#' @param phenotype.type Phenotype covariance structure ('Homo' for homogeneous, 'Heter' for heterogeneous) (default 'Heter').
#' @param lambda Regularization parameter. If NULL, cross-validation will be performed.
#' @param lambda.list List of lambda values to try in cross-validation (used only if lambda is NULL).
#' @param weight Optional weight vector for variance components.
#' @param fixed.effects Optional fixed effects design matrix.
#' @param constrain Logical indicating whether to constrain variance components to be non-negative (default FALSE).
#' @param cv.criteria Criterion for selecting optimal lambda in CV ('LOSS', 'MSE', 'COR') (default is 'LOSS').
#' @param variance.component.list Optional precomputed list of variance components.
#' @param thresh Convergence threshold for iterative methods (default 1e-5).
#' @param max.iter Maximum number of iterations (default 100).
#'
#' @return A list containing estimated variance components (`theta`), phenotype covariance matrix (`phenotype.covariance`),
#'         fixed effects (`beta`), and other matrices.
#' @export
MVKNN = function(phenotype, genotype, input.kernel = NULL, MINQUE.type='MINQUE0', KNN.type='KNN',
                 phenotype.type='Heter', lambda=NULL, lambda.list=NULL, weight=NULL,
                 fixed.effects=NULL, constrain=FALSE, cv.criteria='COR',
                 variance.component.list=NULL, thresh=1e-5, max.iter=100) {

  # If lambda is not specified, perform cross-validation
  if (is.null(lambda)) {
    # Initialize lambda list if not provided
    if (is.null(lambda.list)) {
      lambda.list = c(0.05, 5, 500, 50000, 5000000)
    }

    # Perform cross-validation for MVKNN
    lambda.opt = perform.mvknn.cv(phenotype=phenotype, genotype=genotype, input.kernel=input.kernel,
                                  MINQUE.type=MINQUE.type, KNN.type=KNN.type,
                                  phenotype.type=phenotype.type, lambda.list=lambda.list,
                                  weight=weight, fixed.effects=fixed.effects, constrain=constrain,
                                  criteria=cv.criteria, thresh=thresh, max.iter=max.iter)

    lambda = lambda.opt
  } else {
    lambda.opt = NULL
  }

  # If variance components are not provided, compute them
  if (is.null(variance.component.list)) {
    variance.component.list = KNN2LMM(genotype=genotype, input.kernel=input.kernel, KNN.type=KNN.type)
  }
  n.variance.component = length(variance.component.list)
  sample.size = nrow(phenotype)
  n.traits = ncol(phenotype)
  identity.sample = diag(sample.size)
  identity.trait = diag(n.traits)
  ones.matrix.trait = matrix(1, n.traits, n.traits)

  # Determine phenotype covariance structure
  if (phenotype.type == 'Homo') {
    variance.phenotype = ones.matrix.trait
  } else if (phenotype.type == 'Heter') {
    variance.phenotype = identity.trait
  } else {
    stop("Unsupported phenotype.type.")
  }

  # Initialize fixed effects and residuals
  if (is.null(fixed.effects)) {
    phenotype.vectorized = as.vector(phenotype)
    beta = NULL
  } else {
    n.fixed.effects = ncol(fixed.effects)
    fixed.effects.kronecker = kronecker(identity.trait, fixed.effects)
    phenotype.vectorized = as.vector(phenotype)
    linear.model = lm(phenotype.vectorized ~ fixed.effects.kronecker - 1)
    beta = matrix(coef(linear.model), nrow=n.fixed.effects, ncol=n.traits)
    phenotype.vectorized = residuals(linear.model)
  }

  # Estimate phenotype covariance matrix
  phenotype.covariance = matrix(NA, n.traits, n.traits)

  # Estimate diagonal elements (variances)
  for (trait.index in 1:n.traits) {
    phenotype.indices = ((trait.index - 1) * sample.size + 1):(trait.index * sample.size)
    result = KNN(phenotype=phenotype.vectorized[phenotype.indices], genotype=genotype,
                 input.kernel=input.kernel, MINQUE.type=MINQUE.type, lambda=lambda,
                 weight=weight, constrain=TRUE, variance.component.list=variance.component.list)
    phenotype.covariance[trait.index, trait.index] = result$theta[1]
  }

  # Estimate off-diagonal elements (covariances)
  for (first.trait in 1:(n.traits - 1)) {
    for (second.trait in (first.trait + 1):n.traits) {
      phenotype.first.indices = ((first.trait - 1) * sample.size + 1):(first.trait * sample.size)
      phenotype.second.indices = ((second.trait - 1) * sample.size + 1):(second.trait * sample.size)
      result = KNN.cross(phenotype.vectorized[phenotype.first.indices],
                         phenotype.vectorized[phenotype.second.indices],
                         input.kernel=input.kernel, variance.component.list=list(identity.sample),
                         MINQUE.type=MINQUE.type, lambda=lambda, weight=weight, constrain=TRUE)
      phenotype.covariance[first.trait, second.trait] = result$theta[1]
      phenotype.covariance[second.trait, first.trait] = phenotype.covariance[first.trait, second.trait]
    }
  }

  # Adjust variance components for multivariate case
  variance.component.list.original = variance.component.list
  variance.component.list[[1]] = kronecker(phenotype.covariance, identity.sample)
  if (n.variance.component > 1){
    for (component.index in 2:n.variance.component) {
      variance.component.list[[component.index]] = kronecker(variance.phenotype, variance.component.list[[component.index]])
    }
  }

  # Initialize theta
  if (is.null(weight)) {
    theta = rep(1 / n.variance.component, n.variance.component)
  } else {
    theta = weight
  }

  # Iterative estimation loop
  for (iteration in 1:max.iter) {
    theta.previous = theta

    # Compute V matrix
    variance.matrix = Reduce(`+`, mapply(`*`, variance.component.list, theta.previous, SIMPLIFY=FALSE))
    npd.result = nearPD(variance.matrix)
    variance.matrix = as.matrix(npd.result$mat)
    precision.matrix = solve(variance.matrix)

    # Compute R*variance and R*phenotype
    rotated.variance.list = list()
    for (component.index in seq_len(n.variance.component)) {
      rotated.variance.list[[component.index]] = precision.matrix %*% variance.component.list[[component.index]]
    }
    rotated.phenotype = precision.matrix %*% phenotype.vectorized

    # Compute U matrix
    u.matrix = numeric(n.variance.component)
    for (component.index in seq_len(n.variance.component)) {
      u.matrix[component.index] = sum(rotated.phenotype *
                                        (variance.component.list[[component.index]] %*% rotated.phenotype))
    }

    # Compute C matrix
    covariance.matrix = matrix(0, n.variance.component, n.variance.component)
    for (row.index in seq_len(n.variance.component)) {
      for (col.index in seq_len(n.variance.component)) {
        covariance.matrix[row.index, col.index] = sum(rotated.variance.list[[row.index]] *
                                                        rotated.variance.list[[col.index]])
      }
    }

    # Apply penalty to C matrix
    penalty.matrix = diag(lambda, n.variance.component)
    penalty.matrix[1, 1] = 0  # No penalty on the first variance component
    covariance.matrix = covariance.matrix + penalty.matrix
    npd.result = nearPD(covariance.matrix)
    covariance.matrix = as.matrix(npd.result$mat)

    # Update theta
    covariance.matrix = solve(covariance.matrix)
    theta = covariance.matrix %*% u.matrix

    # Check convergence
    theta.check = theta
    theta.check[theta.check < 0] = 0
    if (sum(abs(theta.check - theta.previous)) < thresh) break

    if (constrain) theta[theta < 0] = 0
  }

  theta = as.numeric(theta)
  result = list(theta=theta, phenotype.covariance=phenotype.covariance, beta=beta,
                covariance.matrix=covariance.matrix, precision.matrix=precision.matrix,
                variance.component.list=variance.component.list,
                variance.component.list.original=variance.component.list.original)

  if (!is.null(lambda.opt)) {
    result$lambda.opt = lambda.opt
  }

  return(result)
}

############################################################
# MVKNN Cross-Validation Function
############################################################

#' @title MVKNN Cross-Validation Function
#'
#' @description
#' Performs cross-validation for lambda selection in MVKNN models.
#' @export
perform.mvknn.cv = function(phenotype, genotype, input.kernel, MINQUE.type, KNN.type,
                            phenotype.type, lambda.list, weight, fixed.effects, constrain,
                            criteria, thresh, max.iter) {

  sample.size = nrow(phenotype)
  n.traits = ncol(phenotype)
  n.lambda.list = length(lambda.list)

  # Split data for cross-validation
  batch.size = sample.size %/% n.lambda.list
  n.genotype.regions = length(genotype)
  randomized.indices = sample(1:sample.size)
  cv.fold.indices = list()

  for (fold.index in 0:(n.lambda.list - 2)) {
    fold.samples = randomized.indices[(fold.index * batch.size + 1):((fold.index + 1) * batch.size)]
    cv.fold.indices = c(cv.fold.indices, list(fold.samples))
  }
  last.fold.samples = randomized.indices[-c(1:((n.lambda.list - 1) * batch.size))]
  cv.fold.indices = c(cv.fold.indices, list(last.fold.samples))

  # Initialize evaluation metrics
  mse.values = numeric(n.lambda.list)
  correlation.values = numeric(n.lambda.list)
  loss.values = numeric(n.lambda.list)

  # Cross-validation loop
  for (lambda.index in 1:n.lambda.list) {
    train.indices = unlist(cv.fold.indices[-lambda.index])
    test.index = unlist(cv.fold.indices[lambda.index])
    phenotype.train = phenotype[train.indices, , drop=FALSE]
    phenotype.test = phenotype[test.index, , drop=FALSE]

    if (!is.null(fixed.effects)) {
      fixed.effects.train = fixed.effects[train.indices, , drop=FALSE]
      fixed.effects.test = fixed.effects[test.index, , drop=FALSE]
    } else {
      fixed.effects.train = NULL
      fixed.effects.test = NULL
    }

    genotype.train = list()
    genotype.test = list()

    for (region.index in 1:n.genotype.regions) {
      genotype.train.region = genotype[[region.index]][train.indices, , drop=FALSE]
      genotype.test.region = genotype[[region.index]][test.index, , drop=FALSE]
      genotype.train = c(genotype.train, list(genotype.train.region))
      genotype.test = c(genotype.test, list(rbind(genotype.train.region, genotype.test.region)))
    }

    # Fit MVKNN model with internal function to avoid recursion
    result = fit.mvknn.internal(phenotype=phenotype.train, genotype=genotype.train,
                                input.kernel=input.kernel, MINQUE.type=MINQUE.type,
                                KNN.type=KNN.type, phenotype.type=phenotype.type,
                                lambda=lambda.list[lambda.index], weight=weight,
                                fixed.effects=fixed.effects.train, constrain=constrain,
                                thresh=thresh, max.iter=max.iter)

    # Predict
    phenotype.predicted = MVKNN.predict(phenotype=phenotype.train, genotype=genotype.test,
                                        theta=result$theta,
                                        phenotype.covariance=result$phenotype.covariance,
                                        input.kernel=input.kernel, KNN.type=KNN.type,
                                        beta=result$beta,
                                        fixed.effects.train=fixed.effects.train,
                                        fixed.effects.test=fixed.effects.test)

    # Calculate evaluation metrics
    # MSE across all traits
    mse.values[lambda.index] = mean((phenotype.test - phenotype.predicted)^2)

    # Average correlation across traits
    trait.correlations = numeric(n.traits)
    for (trait.index in 1:n.traits) {
      trait.correlations[trait.index] = cor(phenotype.test[, trait.index],
                                            phenotype.predicted[, trait.index])
    }
    correlation.values[lambda.index] = mean(trait.correlations, na.rm=TRUE)
  }

  # Select optimal lambda based on criteria
  if (criteria == 'MSE') {
    optimal.index = which.min(mse.values)
  } else if (criteria == 'COR') {
    optimal.index = which.max(correlation.values)
  } else {
    stop("Unsupported criteria for MVKNN.")
  }

  lambda.optimal = lambda.list[optimal.index]
  if (is.na(lambda.optimal)) {
    lambda.optimal = 0
  }

  return(lambda.optimal)
}

############################################################
# Internal MVKNN Fitting Function (to avoid recursion in CV)
############################################################

fit.mvknn.internal = function(phenotype, genotype, input.kernel, MINQUE.type, KNN.type,
                              phenotype.type, lambda, weight, fixed.effects, constrain,
                              thresh, max.iter) {

  # This is essentially the same as MVKNN but without CV logic
  # The code structure follows the main MVKNN function

  # If variance components are not provided, compute them
  variance.component.list = KNN2LMM(genotype=genotype, input.kernel=input.kernel, KNN.type=KNN.type)
  n.variance.component = length(variance.component.list)
  sample.size = nrow(phenotype)
  n.traits = ncol(phenotype)
  identity.sample = diag(sample.size)
  identity.trait = diag(n.traits)
  ones.matrix.trait = matrix(1, n.traits, n.traits)

  # Determine phenotype covariance structure
  if (phenotype.type == 'Homo') {
    variance.phenotype = ones.matrix.trait
  } else if (phenotype.type == 'Heter') {
    variance.phenotype = identity.trait
  } else {
    stop("Unsupported phenotype.type.")
  }

  # Initialize fixed effects and residuals
  if (is.null(fixed.effects)) {
    phenotype.vectorized = as.vector(phenotype)
    beta = NULL
  } else {
    n.fixed.effects = ncol(fixed.effects)
    fixed.effects.kronecker = kronecker(identity.trait, fixed.effects)
    phenotype.vectorized = as.vector(phenotype)
    linear.model = lm(phenotype.vectorized ~ fixed.effects.kronecker - 1)
    beta = matrix(coef(linear.model), nrow=n.fixed.effects, ncol=n.traits)
    phenotype.vectorized = residuals(linear.model)
  }

  # Estimate phenotype covariance matrix
  phenotype.covariance = matrix(NA, n.traits, n.traits)

  # Estimate diagonal elements (variances)
  for (trait.index in 1:n.traits) {
    phenotype.indices = ((trait.index - 1) * sample.size + 1):(trait.index * sample.size)
    result = fit.internal(phenotype=phenotype.vectorized[phenotype.indices], genotype=genotype,
                          input.kernel=input.kernel, MINQUE.type=MINQUE.type, KNN.type=KNN.type,
                          lambda=lambda, weight=weight, fixed.effects=NULL, constrain=TRUE)
    phenotype.covariance[trait.index, trait.index] = result$theta[1]
  }

  # Estimate off-diagonal elements (covariances)
  for (first.trait in 1:(n.traits - 1)) {
    for (second.trait in (first.trait + 1):n.traits) {
      phenotype.first.indices = ((first.trait - 1) * sample.size + 1):(first.trait * sample.size)
      phenotype.second.indices = ((second.trait - 1) * sample.size + 1):(second.trait * sample.size)
      result = KNN.cross(phenotype.vectorized[phenotype.first.indices],
                         phenotype.vectorized[phenotype.second.indices],
                         input.kernel=input.kernel, variance.component.list=list(identity.sample),
                         MINQUE.type=MINQUE.type, lambda=lambda, weight=weight, constrain=TRUE)
      phenotype.covariance[first.trait, second.trait] = result$theta[1]
      phenotype.covariance[second.trait, first.trait] = phenotype.covariance[first.trait, second.trait]
    }
  }

  # Adjust variance components for multivariate case
  variance.component.list.original = variance.component.list
  variance.component.list[[1]] = kronecker(phenotype.covariance, identity.sample)
  if (n.variance.component > 1){
    for (component.index in 2:n.variance.component) {
      variance.component.list[[component.index]] = kronecker(variance.phenotype,
                                                             variance.component.list[[component.index]])
    }
  }

  # Initialize theta
  if (is.null(weight)) {
    theta = rep(1 / n.variance.component, n.variance.component)
  } else {
    theta = weight
  }

  # Iterative estimation loop
  for (iteration in 1:max.iter) {
    theta.previous = theta

    # Compute V matrix
    variance.matrix = Reduce(`+`, mapply(`*`, variance.component.list, theta.previous, SIMPLIFY=FALSE))
    npd.result = nearPD(variance.matrix)
    variance.matrix = as.matrix(npd.result$mat)
    precision.matrix = solve(variance.matrix)

    # Compute R*variance and R*phenotype
    rotated.variance.list = list()
    for (component.index in seq_len(n.variance.component)) {
      rotated.variance.list[[component.index]] = precision.matrix %*% variance.component.list[[component.index]]
    }
    rotated.phenotype = precision.matrix %*% phenotype.vectorized

    # Compute U matrix
    u.matrix = numeric(n.variance.component)
    for (component.index in seq_len(n.variance.component)) {
      u.matrix[component.index] = sum(rotated.phenotype *
                                        (variance.component.list[[component.index]] %*% rotated.phenotype))
    }

    # Compute C matrix
    covariance.matrix = matrix(0, n.variance.component, n.variance.component)
    for (row.index in seq_len(n.variance.component)) {
      for (col.index in seq_len(n.variance.component)) {
        covariance.matrix[row.index, col.index] = sum(rotated.variance.list[[row.index]] *
                                                        rotated.variance.list[[col.index]])
      }
    }

    # Apply penalty to C matrix
    penalty.matrix = diag(lambda, n.variance.component)
    penalty.matrix[1, 1] = 0  # No penalty on the first variance component
    covariance.matrix = covariance.matrix + penalty.matrix
    npd.result = nearPD(covariance.matrix)
    covariance.matrix = as.matrix(npd.result$mat)

    # Update theta
    covariance.matrix = solve(covariance.matrix)
    theta = covariance.matrix %*% u.matrix

    # Check convergence
    theta.check = theta
    theta.check[theta.check < 0] = 0
    if (sum(abs(theta.check - theta.previous)) < thresh) break

    if (constrain) theta[theta < 0] = 0
  }

  theta = as.numeric(theta)

  return(list(theta=theta, phenotype.covariance=phenotype.covariance, beta=beta))
}

############################################################
# MVKNN Mixture Chi-Square Individual Test
############################################################

#' @title MVKNN Mixture Chi-Square Individual Test
#'
#' @description
#' Performs an individual test for variance components in a multivariate KNN model using a mixture of chi-square distributions.
#'
#' @param phenotype Matrix of phenotypes (N x d), where N is the number of samples and d is the number of traits.
#' @param genotype List of N x p genotype matrices.
#' @param test.index Index of the variance component(s) to test.
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param MINQUE.type Type of MINQUE estimator ('MINQUE0' or 'MINQUE1') (default 'MINQUE0').
#' @param KNN.type Type of KNN model ('KNN', 'LMM') (default 'KNN').
#' @param phenotype.type Phenotype covariance structure ('Homo' or 'Heter') (default 'Heter').
#' @param lambda Regularization parameter (default 0).
#' @param weight Optional weight vector for variance components.
#' @param fixed.effects Optional fixed effects design matrix.
#' @param constrain Logical indicating whether to constrain variance components to be non-negative (default FALSE).
#' @param thresh Convergence threshold for iterative methods (default 1e-5).
#' @param max.iter Maximum number of iterations (default 100).
#'
#' @return A list containing the p-value of the test.
#' @export
MVKNN.ind.test = function(phenotype, genotype, test.index, input.kernel = NULL, MINQUE.type='MINQUE0',
                          KNN.type='KNN', phenotype.type='Heter', lambda=0, weight=NULL,
                          fixed.effects=NULL, constrain=TRUE, thresh=1e-5, max.iter=100) {

  # Estimate variance components under the full model
  result.full = MVKNN(phenotype=phenotype, genotype=genotype, input.kernel=input.kernel,
                      MINQUE.type=MINQUE.type, KNN.type=KNN.type, phenotype.type=phenotype.type,
                      lambda=lambda, weight=weight, fixed.effects=fixed.effects,
                      constrain=constrain, thresh=thresh, max.iter=max.iter)
  covariance.matrix = result.full$covariance.matrix
  theta.full = result.full$theta
  precision.matrix = result.full$precision.matrix
  variance.component.list.full.original = result.full$variance.component.list.original

  # Remove the variance component(s) to test to get the reduced model
  variance.component.list.reduced.original = variance.component.list.full.original[-test.index]
  result.reduced = MVKNN(phenotype=phenotype, genotype=genotype, input.kernel=input.kernel,
                         MINQUE.type=MINQUE.type, KNN.type=KNN.type, phenotype.type=phenotype.type,
                         lambda=lambda, weight=weight, fixed.effects=fixed.effects,
                         constrain=constrain, variance.component.list=variance.component.list.reduced.original,
                         thresh=thresh, max.iter=max.iter)
  theta.reduced = result.reduced$theta

  # Compute sigma H1 and sigma H2
  variance.component.list.full = result.full$variance.component.list
  variance.component.list.reduced = result.reduced$variance.component.list
  sigma.hypothesis.one = Reduce(`+`, mapply(`*`, variance.component.list.reduced, theta.reduced, SIMPLIFY=FALSE))
  sigma.hypothesis.two = Reduce(`+`, mapply(`*`, variance.component.list.full, covariance.matrix[test.index, ], SIMPLIFY=FALSE))

  # Compute sigma H
  sigma.h1.svd = svd(sigma.hypothesis.one)
  sigma.h1.sqrt = sigma.h1.svd$u %*% diag(sqrt(sigma.h1.svd$d)) %*% t(sigma.h1.svd$v)
  sigma.hypothesis = sigma.h1.sqrt %*% precision.matrix %*% sigma.hypothesis.two %*%
    precision.matrix %*% sigma.h1.sqrt

  # Compute eigenvalues and p-value using liu method
  eigenvalues = eigen(sigma.hypothesis, symmetric=TRUE, only.values=TRUE)$values

  # Remove very small eigenvalues for numerical stability
  eigenvalues = eigenvalues[eigenvalues > 1e-15]

  p.value = tryCatch({
    davies(q=theta.full[test.index], lambda=eigenvalues)$Qq
  }, error = function(e) {
    # If Liu fails, use Davies
    liu(q=theta.full[test.index], lambda=eigenvalues)
  })

  return(list(p.value=p.value))
}

############################################################
# MVKNN Mixture Chi-Square Overall Test
############################################################

#' @title MVKNN Mixture Chi-Square Overall Test
#'
#' @description
#' Performs an overall test for variance components in a multivariate KNN model using a mixture of chi-square distributions.
#'
#' @param phenotype Matrix of phenotypes (N x d), where N is the number of samples and d is the number of traits.
#' @param genotype List of N x p genotype matrices.
#' @param test.index Indices of the variance components to test.
#' @param input.kernel List specifying the kernel function(s) for each genotype matrix.
#' @param MINQUE.type Type of MINQUE estimator ('MINQUE0' or 'MINQUE1') (default 'MINQUE0').
#' @param KNN.type Type of KNN model ('KNN', 'LMM') (default 'KNN').
#' @param phenotype.type Phenotype covariance structure ('Homo' or 'Heter') (default 'Heter').
#' @param lambda Regularization parameter (default 0).
#' @param weight Optional weight vector for variance components.
#' @param fixed.effects Optional fixed effects design matrix.
#' @param constrain Logical indicating whether to constrain variance components to be non-negative (default FALSE).
#' @param thresh Convergence threshold for iterative methods (default 1e-5).
#' @param max.iter Maximum number of iterations (default 100).
#'
#' @return A list containing the p-value of the test, test mean, and test weights.
#' @export
MVKNN.overall.test = function(phenotype, genotype, test.index=NULL, input.kernel = NULL,
                              MINQUE.type='MINQUE0', KNN.type='KNN', phenotype.type='Heter',
                              lambda=0, weight=NULL, fixed.effects=NULL, constrain=TRUE,
                              thresh=1e-5, max.iter=100) {

  # Estimate variance components under the full model
  result.full = MVKNN(phenotype=phenotype, genotype=genotype, input.kernel=input.kernel,
                      MINQUE.type=MINQUE.type, KNN.type=KNN.type, phenotype.type=phenotype.type,
                      lambda=lambda, weight=weight, fixed.effects=fixed.effects,
                      constrain=constrain, thresh=thresh, max.iter=max.iter)
  covariance.matrix.full = result.full$covariance.matrix
  theta.full = result.full$theta
  precision.matrix.full = result.full$precision.matrix
  variance.component.list.full.original = result.full$variance.component.list.original

  # Remove the variance component(s) to test to get the reduced model
  if (is.null(test.index)) {
    test.index = 1:length(variance.component.list.full.original)
    test.index = test.index[-1]
  }
  variance.component.list.reduced.original = variance.component.list.full.original[-test.index]
  result.reduced = MVKNN(phenotype=phenotype, genotype=genotype, input.kernel=input.kernel,
                         MINQUE.type=MINQUE.type, KNN.type=KNN.type, phenotype.type=phenotype.type,
                         lambda=lambda, weight=weight, fixed.effects=fixed.effects,
                         constrain=constrain, variance.component.list=variance.component.list.reduced.original,
                         thresh=thresh, max.iter=max.iter)
  theta.reduced = result.reduced$theta

  # Compute sigma H1 and sigma H2
  sigma.hypothesis.two.list = list()
  variance.component.list.full = result.full$variance.component.list
  variance.component.list.reduced = result.reduced$variance.component.list
  sigma.hypothesis.one = Reduce(`+`, mapply(`*`, variance.component.list.reduced, theta.reduced, SIMPLIFY=FALSE))
  test.weight = rep(NA, length(test.index))
  index.counter = 1
  test.mean = 0

  for (component.index in test.index) {
    sigma.hypothesis.two = Reduce(f = "+", x = mapply(FUN = "*", variance.component.list.full,
                                                      covariance.matrix.full[component.index,], SIMPLIFY = FALSE))
    sigma.hypothesis.two.list[[index.counter]] = sigma.hypothesis.two
    test.weight[index.counter] = covariance.matrix.full[component.index, component.index]
    test.mean = test.mean + covariance.matrix.full[component.index, component.index] * theta.full[component.index]
    index.counter = index.counter + 1
  }
  sigma.hypothesis.two = Reduce(f = "+", x = mapply(FUN = "*", sigma.hypothesis.two.list, test.weight, SIMPLIFY = FALSE))

  # Compute sigma H
  sigma.h1.svd = svd(sigma.hypothesis.one)
  sigma.h1.sqrt = sigma.h1.svd$u %*% diag(sqrt(sigma.h1.svd$d)) %*% t(sigma.h1.svd$v)
  sigma.hypothesis = sigma.h1.sqrt %*% precision.matrix.full %*% sigma.hypothesis.two %*%
    precision.matrix.full %*% sigma.h1.sqrt

  # Compute eigenvalues and p-value using liu method
  eigenvalues = eigen(sigma.hypothesis, symmetric=TRUE, only.values=TRUE)$values

  # Remove very small eigenvalues for numerical stability
  eigenvalues = eigenvalues[eigenvalues > 1e-15]

  p.value = tryCatch({
    davies(q = test.mean, lambda = eigenvalues)$Qq
  }, error = function(e) {
    # If Liu fails, use Davies
    liu(q = test.mean, lambda = eigenvalues)
  })

  return(list(p.value=p.value, test.mean=test.mean, test.weight=test.weight))
}
