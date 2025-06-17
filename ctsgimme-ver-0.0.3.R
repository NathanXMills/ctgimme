#nolint start
# Why are do so many of these have defualt parameters of Null? It seems like this would break things later?
ctsgimme = function(varnames = NULL, dataframe = NULL,
                    id = NULL, time = NULL,
                    cores = NULL, directory = NULL, 
                    sig.thrsh = 0.75,
                    ben.hoch = TRUE,
                    Galpha = 0.05,
                    Ialpha = 0.01,
                    ME.var = diag(1e-5, length(varnames)), 
                    PE.var = diag(1.00, length(varnames)),
                    time.intervals = c(1)){
  # figures out which connection should be the next one made; which one would best improve the model
  JPmx = function (model, matrices = NA, full = TRUE){
    OpenMx:::warnModelCreatedByOldVersion(model)
    # Checks to see if the matrices are NA (you did not provide one), and if so then it creates them from the model # nolint: line_length_linter.
    if (OpenMx:::single.na(matrices)) {
      matrices = names(model$matrices)
      #if the model has things that it wants to be exempted from the matrices this will remove them  # nolint: line_length_linter.
      #Unsure what MXExpectationRAM is or why model$exctation would be that, unless it's just something that you should know and is set by person making the model # nolint: line_length_linter.
      if (is(model$expectation, "MxExpectationRAM")) {
        matrices = setdiff(matrices, model$expectation$F)
      }
    }
    if (imxHasWLS(model)) {
      stop("modification indices not implemented for WLS fitfunction")
    }
    param = omxGetParameters(model)
    param.names = names(param)
      # makes a copy of the model but has every parameter fixed
    gmodel = omxSetParameters(model, free = FALSE, labels = param.names)
    mi.r = NULL
    mi.f = NULL
    EPCs = NULL
    a.names = NULL
    new.models = list()
    for (amat in matrices) {
      # does this break if you directly give matrices that are not inline with the model?
      # what object type is model? What does it contain aside from matrices?
      matObj = model[[amat]]
        # This is the same demensions as the matrix, but each cell is true/false
        # Used later to keep track of which cells have been tried freed
      freemat = matObj$free
        # unsure what this half empty matrix is for?
      sym.sel = upper.tri(freemat, diag = TRUE)
      notSymDiag = !(is(gmodel[[amat]])[1] %in% c("DiagMatrix", 
                                                  "SymmMatrix"))
      for (i in 1:length(freemat)) {
          # checks to see if we have already tried freeing this cell
        if (freemat[i] == FALSE && (notSymDiag || sym.sel[i] == 
                                    TRUE)) {
          tmpLab = gmodel[[amat]]$labels[i]
          plusOneParamModel = model
            # frees the cell that we are looking at if it's not empty or null
          if (length(tmpLab) > 0 && !is.na(tmpLab)) {
            gmodel = omxSetParameters(gmodel, labels = tmpLab, 
                                      free = TRUE)
            plusOneParamModel = omxSetParameters(plusOneParamModel, 
                                                 labels = tmpLab, free = TRUE)
          }
            # frees the cell manually
            # why does this distinction exist? Why can't you just use the else here?
          else {
            gmodel[[amat]]$free[i] = TRUE
            plusOneParamModel[[amat]]$free[i] = TRUE
          }
          # I don't understand why Zero, Diagnal, and Symmetrical matricies don't work? Also why symmetrical -> symmetrical?
            # turns zero matrices into full matricies, and diagnal and symmetrical ones into symmetrical matricies
          if (is(gmodel[[amat]])[1] %in% c("ZeroMatrix")) {
            cop = gmodel[[amat]]
            newSingleParamMat = mxMatrix("Full", nrow = nrow(cop), 
                                         ncol = ncol(cop), values = cop$values, free = cop$free, 
                                         labels = cop$labels, name = cop$name, lbound = cop$lbound, 
                                         ubound = cop$ubound, dimnames = dimnames(cop))
            bop = plusOneParamModel[[amat]]
            newPlusOneParamMat = mxMatrix("Full", nrow = nrow(bop), 
                                          ncol = ncol(bop), values = bop$values, free = bop$free, 
                                          labels = bop$labels, name = bop$name, lbound = bop$lbound, 
                                          ubound = bop$ubound, dimnames = dimnames(bop))
          }
          else if (is(gmodel[[amat]])[1] %in% c("DiagMatrix", 
                                                "SymmMatrix")) {
            cop = gmodel[[amat]]
            newSingleParamMat = mxMatrix("Symm", nrow = nrow(cop), 
                                         ncol = ncol(cop), values = cop$values, free = (cop$free | 
                                                                                          t(cop$free)), labels = cop$labels, name = cop$name, 
                                         lbound = cop$lbound, ubound = cop$ubound, 
                                         dimnames = dimnames(cop))
            bop = plusOneParamModel[[amat]]
            newPlusOneParamMat = mxMatrix("Symm", nrow = nrow(bop), 
                                          ncol = ncol(bop), values = bop$values, free = (bop$free | 
                                                                                           t(bop$free)), labels = bop$labels, name = bop$name, 
                                          lbound = bop$lbound, ubound = bop$ubound, 
                                          dimnames = dimnames(bop))
          }
          else {
            newSingleParamMat = gmodel[[amat]]
            newPlusOneParamMat = plusOneParamModel[[amat]]
          }
          gmodel[[amat]] = newSingleParamMat
          plusOneParamModel[[amat]] = newPlusOneParamMat
            # I don't understand what these next couple lines are doing. 
            # I think they're finding some things that are used for the modification index, but i don't know what those are
          custom.compute = mxComputeSequence(list(mxComputeNumericDeriv(checkGradient = FALSE), 
                                                  mxComputeReportDeriv()))
          gmodel = mxModel(gmodel, custom.compute)
          grun = try(mxRun(gmodel, silent = TRUE, suppressWarnings = FALSE, 
                           unsafe = TRUE))
          nings = TRUE
          if (is(grun, "try-error")) {
            gmodel = omxSetParameters(gmodel, labels = names(omxGetParameters(gmodel)), 
                                      free = FALSE)
            next
          }
          # this calculates the modification index which tells us if changing this cell would improve the model 
          grad = grun$output$gradient
          hess = grun$output$hessian
          modind = 0.5 * grad^2/hess
            # this calculates a more accurate modification index if you want it
            # why not always use this one? why even have the above calculation?
          if (full == TRUE) {
            custom.compute.smart = mxComputeSequence(list(mxComputeNumericDeriv(knownHessian = model$output$hessian, 
                                                                                checkGradient = FALSE), mxComputeReportDeriv()))
            plusOneParamRun = mxRun(mxModel(plusOneParamModel, 
                                            custom.compute.smart), silent = TRUE, suppressWarnings = FALSE, 
                                    unsafe = TRUE)
            grad.full = plusOneParamRun$output$gradient
            grad.full[is.na(grad.full)] = 0
            hess.full = plusOneParamRun$output$hessian
            modind.full = 0.5 * t(matrix(grad.full)) %*%
              solve(hess.full) %*% matrix(grad.full)
            EPC = modind.full/grad.full[grad.full!=0]
          }
          else {
            modind.full = NULL
          }
          # if this cell should be changed it records it
          n.names = names(omxGetParameters(grun))
          if (length(modind) > 0) {
            a.names = c(a.names, n.names)
            mi.r = c(mi.r, modind)
            mi.f = c(mi.f, modind.full)
            EPCs = c(EPCs, EPC)
            new.models = c(new.models, plusOneParamModel)
          }
          # unfrees the current parameter before the next one is tested
          gmodel = omxSetParameters(gmodel, labels = names(omxGetParameters(gmodel)), 
                                    free = FALSE)
        }
      }
      names(mi.r) = a.names
      if (full == TRUE) {
        names(mi.f) = a.names
        names(EPCs) = a.names
      }
      names(new.models) = a.names
    }
    # repeats the process for any and all submodels
      # I don't know what submodels these are. Are they subgroups, or individuals inside of group level analysis? 
    if (length(model$submodels) > 0) {
      for (asubmodel in names(model$submodels)) {
        ret = c(ret, JPmx(asubmodel))
      }
    }
    # message(paste0("The N is ", N))
    return(list(MI = mi.r, MI.Full = mi.f, plusOneParamModels = new.models, EPC = -EPCs))
  }
  # Fixed to 1.00 until official subgrouping is tested.
  sub.sig.thrsh = 1.00
  S.Galpha = 0.05
  #  Creating Directory---###---###---###
  ###---###---###---###---###---###---###
  dir.create(directory, showWarnings = FALSE)
  dir.create(paste0(directory, "/MIs/"), showWarnings = FALSE)
  dir.create(paste0(directory, "/Models/"), showWarnings = FALSE)
  dir.create(paste0(directory, "/Models/Individuals/"), showWarnings = FALSE)  
    # What does this dataframe look like? Trouble looking at Dataset.RDS 
  # Specifying Cores
  ids = unique(dataframe[,id])
  nvar = length(varnames)
    # In theory shouldn't it be possible for there to be no group level commonalities? (maybe incase of null?)
    # Why do you need to know the number of cores before, and why should they be based on the numher of ids?
    # Shouldn't GIMME find these automatically, why are we getting this number from number of ids?
  if(length(ids) < cores){
    cores = length(ids)
    message("Cores adjusted to your sample-size!")
  }else{cores = cores}
  ###---###---###---###---###---###---###
  # Setting Up Multiple Model Estimation#
  ###---###---###---###---###---###---###
  library(parallel)
  library(dynr)
  library(OpenMx)
  library(igraph)
  library(qgraph)
    # I believe this is setting up parrelel computing for individual analysis
  cl = makeCluster(cores, type = "PSOCK")
  clusterExport(cl, c("dataframe", "JPmx",
                      "varnames", "id", "time", "directory",
                      "nvar", "ME.var", "PE.var"),
                envir = environment())
  clusterEvalQ(cl, {
    packs = list('dynr', 'OpenMx', 'qgraph')
    invisible(lapply(packs, require, character.only = T))
  })
  ###---###---###---###---###---###---###
  # Running Step 1 in Parallel---###---##
  ###---###---###---###---###---###---###
    # this looks at individuals before the group level analysis, and creates the models for what everperson looks like
    # i is whatever id is currently being looked at, simular to looping through a list. Just it's all simultanious
  parLapply(cl, ids, function(i) {
      # Creates a matrixs with the diagnol being "A_1,1", "A_2,2" etc
    DRIFT = diag(paste0("A_", 1:nvar, ",", 1:nvar), nvar)
    # just looking at the part of the dataframe related to the current id
    subset_dat = subset(dataframe, id == i)
    # subset_dat[,varnames] = scale(subset_dat[,varnames])
      # I don't realy understand what all of these do
    amat = mxMatrix("Full", nvar, nvar,
                    free   = DRIFT != "0",
                    name   = "A", ubound = 30, lbound = -30)
    bmat = mxMatrix('Zero', nvar, nvar, name='B')
    cdim = list(varnames, paste0('F', 1:nvar))
    cmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name = 'C', dimnames = cdim)
    dmat = mxMatrix('Zero', nvar, nvar, name='D')
    qmat = mxMatrix('Symm', nvar, nvar, FALSE, PE.var, name='Q')
    rmat = mxMatrix('Symm', nvar, nvar, FALSE, ME.var, name='R')
    xmat = mxMatrix('Full', nvar, 1, FALSE, rep(0, nvar), name='x0')
    pmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name='P0')
    umat = mxMatrix('Zero', nvar, 1, name='u')
    tmat = mxMatrix('Full', 1, 1, FALSE, name='time', labels='data.Time')
      # creates a ct OMX model based on the privious data
    osc = mxModel("OUMod", 
                  amat, bmat, cmat, dmat, qmat, 
                  rmat, xmat, pmat, umat, tmat,
                  mxExpectationSSCT('A', 'B', 'C', 'D', 'Q', 
                                    'R', 'x0', 'P0', 'u', 'time'),
                  mxFitFunctionML(),
                  mxData(subset_dat, 'raw'))
      # fits the model using mxTryHard()
    analysis_result = tryCatch({
      fit = mxTryHard(osc)
    }, error = function(e) {
      message("Error for subject ", i, ": ", e$message)
    })
      # figures out what is the next connection that should be made to improve the model and makes it
    if (!is.null(analysis_result)) {
      MIs = JPmx(analysis_result, matrices = "A")
      saveRDS(
        object = MIs,
        file   = paste0(directory, "/MIs/MI_", i, ".RDS")
      )
    }
  })
    # stops the parralel processing now that individual analysis is done
  stopCluster(cl)
  ###---###---###---###---###---###---###
  # Iteration Steps
  ###---###---###---###---###---###---###
  iterate = 0; count = 1
  m = (nvar^2)-nvar
  ks = matrix(NA, m, 1)
  for(k in 1:m){
    ks[k,] = (k/m) * Galpha
  }
  ks = matrix(sort(ks, TRUE), m, 1)
  if(ben.hoch == FALSE){
    ks = matrix(Galpha, nrow(ks), 1)
  }
  DRIFT = diag(paste0("A_", 1:nvar, 1:nvar), nvar)
    # this is the core of the group level grouping and connections
  while(iterate < 1){
      # reads all of the previously made modification indices and stores them in "files"
    rdss = list.files(paste0(directory, "/MIs/"), pattern = "\\.RDS$", full.names = TRUE)
    files = NULL
    for (file in rdss) {
      file_id = gsub("MI_|\\.RDS", "", basename(file))
      files = cbind(files, tryCatch({
        c(readRDS(file)$"MI.Full")
      }, error = function(e) {
        message("Failed to read ", file, ": ", e$message)
        NULL
      })) 
    }
      # finds the best, most significant, possible change. Then sees if it is a significant change and checks if it is above the set threshold
    SigThresh = sum(pchisq(files[which.max(rowMeans(files)),], 1, lower.tail = FALSE) < ks[count,]) / 
      ncol(files) >= sig.thrsh
    
      # adds the connection to the group level if it met the requirements
      # deletes all modification indices
    if(SigThresh){
      param.to.add = which.max(rowMeans(files))
      cells = as.numeric(unlist(regmatches(rownames(files)[param.to.add], 
                                           gregexpr("\\d+", rownames(files)[param.to.add])))) # nolint
      DRIFT[cells[1], cells[2]] = paste0("A_", cells[1], ",", cells[2])
      message(paste0("Adding drift parameter A[", cells[1], ",", cells[2],"]"))
      message(paste0("Completed Step ", count))
        # I don't understand, why do we do this 4 time?
      unlink(paste0(directory, "/MIs/", "*"), recursive = FALSE, force = TRUE)
      unlink(paste0(directory, "/MIs/", "*"), recursive = FALSE, force = TRUE)
      unlink(paste0(directory, "/MIs/", "*"), recursive = FALSE, force = TRUE)
      unlink(paste0(directory, "/MIs/", "*"), recursive = FALSE, force = TRUE)
      count = count + 1
    }else{
        # Saves the current group level model and leaves the loop
      G.DRIFT = DRIFT
      output_path = file.path(directory, "Group Paths.png")
      png(filename = output_path, width = 800, height = 800)
      # Plot the community structure
      
      qgraph(t((G.DRIFT != "0") * 1), layout = "circle", labels = varnames, 
             edge.width = 5, diag = TRUE, edge.labels = "GROUP")
      dev.off()
      ks = ks[c(count:nrow(ks)),]
      message("Group Search Complete.")
      iterate = 1
    }
      # Reruns and remakes individual models based on the new group level models
      # simular to lines 211-252, possible to condense them?
    cl = makeCluster(cores, type = "PSOCK")
    clusterExport(cl, c("dataframe", "JPmx", "ME.var", "PE.var",
                        "varnames", "id", "time", "directory",
                        "nvar", "DRIFT", "count"),
                  envir = environment())
    clusterEvalQ(cl, {
      packs = list('dynr', 'OpenMx', 'qgraph')
      invisible(lapply(packs, require, character.only = T))
    })
    parLapply(cl, ids, function(i) {
      subset_dat = subset(dataframe, id == i)
      subset_dat[,varnames] = scale(subset_dat[,varnames])
      amat = mxMatrix("Full", nvar, nvar,
                      free   = DRIFT != "0",
                      name   = "A", ubound = 30, lbound = -30)
      bmat = mxMatrix('Zero', nvar, nvar, name='B')
      cdim = list(varnames, paste0('F', 1:nvar))
      cmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name = 'C', dimnames = cdim)
      dmat = mxMatrix('Zero', nvar, nvar, name='D')
      qmat = mxMatrix('Symm', nvar, nvar, FALSE, PE.var, name='Q')
      rmat = mxMatrix('Symm', nvar, nvar, FALSE, ME.var, name='R')
      xmat = mxMatrix('Full', nvar, 1, FALSE, rep(0, nvar), name='x0')
      pmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name='P0')
      umat = mxMatrix('Zero', nvar, 1, name='u')
      tmat = mxMatrix('Full', 1, 1, FALSE, name='time', labels='data.Time')
      osc = mxModel("OUMod", 
                    amat, bmat, cmat, dmat, qmat, 
                    rmat, xmat, pmat, umat, tmat,
                    mxExpectationSSCT('A', 'B', 'C', 'D', 'Q', 
                                      'R', 'x0', 'P0', 'u', 'time'),
                    mxFitFunctionML(),
                    mxData(subset_dat, 'raw'))  
      analysis_result = tryCatch({
        fit = mxTryHard(osc)
      }, error = function(e) {
        message("Error for subject ", i, ": ", e$message)
      })
      if (!is.null(analysis_result)) {
        MIs = JPmx(analysis_result, matrices = "A")
        saveRDS(
          object = MIs,
          file   = paste0(directory, "/MIs/", "/MI_", i, ".RDS")
        )
        saveRDS(
          object = analysis_result,
          file   = paste0(directory, "/Models/", "/Model_", i, ".RDS")
        )
      }
    })
    stopCluster(cl)
  }
  ###---###---###---###---###---###---###
  # Subgrouping Stage
  ###---###---###---###---###---###---###
  if(sub.sig.thrsh == 1.00){
    memb = rep(1, length(ids))
    message("Subgrouping Disabled. Undergoing Testing and Evaluation")
  }
  
  ###---###---###---###---###---###---###
  # Iteration Steps - Subgroup
  ###---###---###---###---###---###---###
  for(subgroup in sort(unique(memb))){
    # matrix showing what parameters are freed
      # I think this came from the group level analysis, but unsure why you would use it here?
      # Maybe it is for telling which cells are common at group level, so not to touch those?
    DRIFT = G.DRIFT
    dir.create(paste0(directory, "/Models/Subgroup ", subgroup, "/"), showWarnings = FALSE)
    # The following section has a ton of commonalities to lines 259+, and also another section later. Possible place to combine?
    iterate = 0; count = 1
    m = (nvar^2)-nvar
    sg.ks = matrix(NA, m, 1)
    for(k in 1:m){
      sg.ks[k,] = (k/m) * S.Galpha
    }
    sg.ks = matrix(sort(sg.ks, TRUE), m, 1)
    if(ben.hoch == FALSE){
      sg.ks = matrix(S.Galpha, nrow(sg.ks), 1)
    }
    memb.id = cbind(unique(dataframe$id), memb)
    while(iterate < 1){
      # this first part is VERY simular to lines 271-296ish. Only real change is how MIs are gotten
      # Finds the data for all individuals in the subgroup, then finds their modification indices
      new.data = subset(dataframe, id %in% subset(memb.id[,1], memb == subgroup))
      rdss = list.files(paste0(directory, "/MIs/"), pattern = "\\.RDS$", full.names = TRUE)
      rds_ids = as.numeric(gsub("MI_|\\.RDS", "", basename(rdss)))
      valid_ids = unique(new.data$id)
      matching_indices = which(rds_ids %in% valid_ids)
      matching_files = rdss[matching_indices][order(rds_ids[matching_indices])]
      files = NULL
      # Combines all the files for the individuals in the subgroup
      for (file in matching_files) {
        file_id = gsub("MI_|\\.RDS", "", basename(file))
        files = cbind(files, tryCatch({
          c(readRDS(file)$"MI.Full")
        }, error = function(e) {
          message("Failed to read ", file, ": ", e$message)
          NULL
        }))
      }
        # finds the best, most significant, possible change. Then sees if it is a significant change and checks if it is above the set threshold
      SigThresh = sum(pchisq(files[which.max(rowMeans(files)),], 1, lower.tail = FALSE) < sg.ks[count,]) / 
        ncol(files) >= sig.thrsh
      
      if(SigThresh){
        param.to.add = which.max(rowMeans(files))
        cells = as.numeric(unlist(regmatches(rownames(files)[param.to.add], 
                                             gregexpr("\\d+", rownames(files)[param.to.add]))))
        DRIFT[cells[1], cells[2]] = paste0("A_", cells[1], ",", cells[2])
        message(paste0("Adding drift parameter A[", cells[1], ",", cells[2],"]"))
        message(paste0("Completed Step ", count))
        count = count + 1
      }else{
        # Finished the subgroup, either move to individual level or do other subgroups
        if(max(memb) == 1){
          message("Beginning Individual-Level Search.")
        }else{
          message(paste0("Subgroup Search ", subgroup," of ", max(memb)," Complete."))
          output_path = file.path(paste0(directory, "/Models/Subgroup ", subgroup, "/Subgroup ", subgroup, " Paths.png"))
          png(filename = output_path, width = 800, height = 800)
          # Plot the subgroup structure
          qgraph(t(abs(((G.DRIFT != "0") * 1) - ((DRIFT != "0") * 1))), layout = "circle", labels = varnames, 
                 edge.width = 5, diag = TRUE, edge.labels = paste0("SG-", subgroup))
          dev.off()
          message(paste0("Beginning Individual Model Fitting for Subgroup Members."))
        }
          # The following section has a ton of commonalities to lines 259+, and also another section
        m = (nvar^2)-sum(DRIFT!="0")
        nks = matrix(NA, m, 1)
        for(k in 1:m){
          nks[k,] = (k/m) * Ialpha
        }
        nks = matrix(sort(nks, TRUE), m, 1)
        SG.DRIFT = DRIFT
        cl = makeCluster(cores, type = "PSOCK")
        clusterExport(cl, c("new.data", "JPmx", "ME.var", "PE.var",
                            "varnames", "id", "time", "directory",
                            "nvar", "SG.DRIFT", "count",
                            "subgroup", "nks"),
                      envir = environment())
        clusterEvalQ(cl, {
          packs = list('dynr', 'OpenMx', 'qgraph')
          invisible(lapply(packs, require, character.only = T))
        })
        parLapply(cl, valid_ids, function(i) {
          DRIFT = SG.DRIFT
          subset_dat = subset(new.data, id == i)
          amat = mxMatrix("Full", nvar, nvar,
                          free   = DRIFT != "0",
                          name   = "A", ubound = 30, lbound = -30)
          bmat = mxMatrix('Zero', nvar, nvar, name='B')
          cdim = list(varnames, paste0('F', 1:nvar))
          cmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name = 'C', dimnames = cdim)
          dmat = mxMatrix('Zero', nvar, nvar, name='D')
          qmat = mxMatrix('Symm', nvar, nvar, FALSE, PE.var, name='Q')
          rmat = mxMatrix('Symm', nvar, nvar, FALSE, ME.var, name='R')
          xmat = mxMatrix('Full', nvar, 1, FALSE, rep(0, nvar), name='x0')
          pmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name='P0')
          umat = mxMatrix('Zero', nvar, 1, name='u')
          tmat = mxMatrix('Full', 1, 1, FALSE, name='time', labels='data.Time')
          osc = mxModel("OUMod", 
                        amat, bmat, cmat, dmat, qmat, 
                        rmat, xmat, pmat, umat, tmat,
                        mxExpectationSSCT('A', 'B', 'C', 'D', 'Q', 
                                          'R', 'x0', 'P0', 'u', 'time'),
                        mxFitFunctionML(),
                        mxData(subset_dat, 'raw'))  
          fit = mxTryHard(osc)
          fit2 = fit
          optimization = 0
          count = 1
          while(optimization < 1){
            if(count > 1){
              fit2 = fit
              fit = mxTryHard(osc)
              if(any(is.na(summary(fit)$parameters$Std.Error))){
                fit = fit2
                optimization = 1
                break
              }
            }
            MIs = JPmx(fit, matrices = "A")
            if(is.null(MIs)){optimization = 1; fit = fit2}
            if(abs(MIs$MI.Full)[which.max(abs(MIs$MI.Full))] >= qchisq(1-nks[count,], df = 1)){
              cells = as.numeric(unlist(regmatches(names(which.max(MIs$MI.Full)),
                                                   gregexpr("\\d+", names(which.max(MIs$MI.Full))))))
              osc$A$free[cells[1], cells[2]] = TRUE
              osc$A$labels[cells[1], cells[2]] = paste0("A_", cells[1], ",", cells[2])
              message(paste0("Adding drift parameter A[", cells[1], ",", cells[2],"]"))
              MIs = NULL
              count = count + 1
              if(sum(osc$A$free) == nvar^2){
                optimization = 1
              }
            }else{
              optimization = 1
            }
          }
          message("Pruning Stage.")
          prune = 0
          count = 1
          stat.sig1 = summary(fit)$parameters
          while (prune < 1) {
            stat.sig = summary(fit)$parameters
            ests = matrix(0, nvar, nvar)
            for(jj in 1:nrow(stat.sig)){
              ests[stat.sig[,3][jj], stat.sig[,4][jj]] = stat.sig[jj,5]
            }
            prunable = stat.sig[stat.sig[,5] == ests[setdiff(which(ests != 0), which(ests != 0 & DRIFT != "0"))],]
            if (nrow(prunable) == 0 | is.null(nrow(prunable))) {
              message("No unprotected parameters left to prune.")
              break
            }
            se = prunable[, 6] * qnorm(0.975)
            z.scores = abs(prunable[, 5]) / se
            min_z_index = which.min(z.scores)
            if(length(z.scores[min_z_index]) == 0){
              message("All Prunable Paths Removed.")
              break
              prune = 1
            }
            if (z.scores[min_z_index] < qnorm(0.975)) {
              this = prunable[min_z_index,]
              cells = matrix(c(this$row, this$col), 1, 2)
              osc$A$free[cells[1,1], cells[1,2]] = FALSE
              fit = mxTryHard(osc)
              message(paste0("NOTE: Pruned drift parameter: A[", cells[1,1], ",", cells[1,2], "]!"))
            } else {
              message("No Pruning Conducted.")
              prune = 1
            }
          }
          if (!is.null(fit)) {
            saveRDS(
              object = fit,
              file   = paste0(directory, "/Models/Individuals/FinalModel_", i, ".RDS")
            )
            
            stat.sig = summary(fit)$parameters
            ests = matrix(0, nvar, nvar)
            for(jj in 1:nrow(stat.sig)){
              ests[stat.sig[,3][jj], stat.sig[,4][jj]] = stat.sig[jj,5]
            }
            png(filename = paste0(directory, "/Models/Individuals/FinalModel_", i, ".PNG"), 
                width = 800, height = 800)
            qgraph(t(ests), layout = "circle", labels = varnames, 
                   edge.width = 1, diag = TRUE, edge.labels = round(c(t(ests)), 2),
                   theme = "colorblind", fade = FALSE)
            dev.off()
            
            unlink(list.files(paste0(directory, "/Models/"), pattern = "\\.RDS$", full.names = TRUE))
            unlink(list.files(paste0(directory, "/Models/"), pattern = "\\.RDS$", full.names = TRUE))
            unlink(list.files(paste0(directory, "/Models/Subgroup ",subgroup,"/"), pattern = "\\.RDS$", full.names = TRUE))
          }
        })
        stopCluster(cl)
        iterate = 1
      }
      cl = makeCluster(cores, type = "PSOCK")
      clusterExport(cl, c("new.data", "JPmx", "ME.var", "PE.var",
                          "varnames", "id", "time", "directory",
                          "nvar", "DRIFT", "count",
                          "subgroup"),
                    envir = environment())
      clusterEvalQ(cl, {
        packs = list('ctsem', 'ctsemOMX', 'dynr', 'OpenMx', 'qgraph')
        invisible(lapply(packs, require, character.only = T))
      })
      qmat = mxMatrix('Diag', nvar, nvar, FALSE, PE.var, name='Q')
      rmat = mxMatrix('Diag', nvar, nvar, FALSE, ME.var, name='R')
      parLapply(cl, valid_ids, function(i) {
        subset_dat = subset(new.data, id == i)
        amat = mxMatrix("Full", nvar, nvar,
                        free   = DRIFT != "0",
                        name   = "A", ubound = 30, lbound = -30)
        bmat = mxMatrix('Zero', nvar, nvar, name='B')
        cdim = list(varnames, paste0('F', 1:nvar))
        cmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name = 'C', dimnames = cdim)
        dmat = mxMatrix('Zero', nvar, nvar, name='D')
        qmat = mxMatrix('Symm', nvar, nvar, FALSE, PE.var, name='Q')
        rmat = mxMatrix('Symm', nvar, nvar, FALSE, ME.var, name='R')
        xmat = mxMatrix('Full', nvar, 1, FALSE, rep(0, nvar), name='x0')
        pmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name='P0')
        umat = mxMatrix('Zero', nvar, 1, name='u')
        tmat = mxMatrix('Full', 1, 1, FALSE, name='time', labels='data.Time')
        osc = mxModel("OUMod", 
                      amat, bmat, cmat, dmat, qmat, 
                      rmat, xmat, pmat, umat, tmat,
                      mxExpectationSSCT('A', 'B', 'C', 'D', 'Q', 
                                        'R', 'x0', 'P0', 'u', 'time'),
                      mxFitFunctionML(),
                      mxData(subset_dat, 'raw'))  
        analysis_result = tryCatch({
          fit = mxTryHard(osc)
        }, error = function(e) {
          message("Error for subject ", i, ": ", e$message)
        })
        if (!is.null(analysis_result)) {
          MIs = JPmx(analysis_result, matrices = "A")
          saveRDS(
            object = MIs,
            file   = paste0(directory, "/MIs/", "/MI_", i, ".RDS")
          )
          saveRDS(
            object = analysis_result,
            file   = paste0(directory, "/Models/Subgroup ", subgroup, "/Model_", i, ".RDS")
          )
        }
      })
      stopCluster(cl)
    }
  }
  message(paste0("Subgrouping with Continuous-Time GIMME Complete. Find networks in ", directory, "."))
  unlink(file.path(directory, "MIs"), recursive = TRUE, force = TRUE)
  unlink(file.path(directory, "MIs"), recursive = TRUE, force = TRUE)
  unlink(file.path(directory, "MIs"), recursive = TRUE, force = TRUE)
  unlink(file.path(directory, "MIs"), recursive = TRUE, force = TRUE)
  if(sub.sig.thrsh == 1){
    return(message("Continuous-Time S-GIMME Complete."))
  }else{return(walktrap_comm)}
}

#nolint end