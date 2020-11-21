#Functions to grab a plan and rename all of its targets +
#external targets used as arguments.


#Rename a single argument
#If the argument is a call, loop modify_call_args. Allows for renaming arguments within nested calls
#@param: arg - argument, will be modified if of class 'name' or 'call'
#in_regex: regex expression to identify target name within argument
#in_suffix: suffix to add to target name

modify_arg <- function(arg, in_regex, in_suffix, verbose = FALSE) {
  if (inherits(arg, 'name')) {
    argc <- as.character(arg) #Convert argument value to character
    argc_match <- grep(in_regex, argc, value = T, perl = T) #Check whether argument value it is a target name
    if (length(argc_match)==1) { #If so
      argrep <- gsub(in_regex, paste0(argc_match, in_suffix), argc, perl=T) #Add suffix to target name within argument value
      if (verbose) {
        print(paste('Changing argument', argc, 'to', argrep)) #Print change
      }
      return(as.name(argrep)) #Return argument value as object of class 'name' (same as input)
    }
  } else if (inherits(arg, 'call')) { #If argument value is of class 'call'
    #Run nested modification of all argument values
    argrep <- modify_call_args(call=arg,
                               in_regex = in_regex,
                               in_suffix = in_suffix)
    return(argrep)
  }
}

#Modify all arguments within a call.
#call: call whose argument values must be modified
#in_regex: regular expression to recognize target names within argument values
#in_suffix: string to add at the end of target name

modify_call_args <- function(call, in_regex, in_suffix, verbose = FALSE) {
  if (inherits(call, 'call')) { #If 'call' actually of class call
    args <- rlang::call_args(call)  #Get call arguments

    args_modified <- lapply(args, modify_arg,
                           in_regex=in_regex, in_suffix=in_suffix, verbose = verbose) #Modify all argument values
    arg_names <- rlang::call_args_names(call) #Get argument names

    if (any(arg_names=="") | any(duplicated(arg_names))) { #If there are duplicate argument names or unnamed arguments (as in a primitive function like foo$bar)
      #Keep the original argument value if not modified
      same_i <- which(unlist(lapply(args_modified, is.null)))
      if (length(same_i) > 0) {
        for (i in same_i) {
          args_modified[[i]] <- args[[i]]
        }
      }

      #Recreate call
      call_modified <- rlang::call2(rlang::call_name(call), !!!args_modified)

    } else { #if all arguments are named and unique
      call_modified <- rlang::call_modify(call, !!!plyr::compact(args_modified))
    }

    return(call_modified)
  } else { #If 'call' is not of class call (e.g. if it is e.g. just a string or number)
    return(call)
  }
}

#Rename all targets of a plan + external targets used as arguments.
#@plan: to modify
#@branch_suffix: character to add at the end of every target name
#@external_arguments_to_modify: character vector of targets used as arguments in plan to modify as well

branch_plan <- function(plan, branch_suffix, external_arguments_to_modify=NULL, verbose = FALSE) {
  plan_modif <- plan #make a copy of plan

  #Create regular expression to recognize all target names (in addition to "external_arguments_to_modify")
  target_regex <- paste(paste0('(',
                               c(plan$target,
                                 external_arguments_to_modify),
                               '(?=([$].+)*))$'),
                        collapse='|')
  #Add suffix to all target names
  plan_modif$target <- paste0(plan_modif$target, branch_suffix)

  #Add suffix to all target names within commands/calls
  plan_modif$command <- lapply(
    plan_modif$command, function(command) {
      #print(command)
      modify_call_args(call=command,
                       in_regex = target_regex,
                       in_suffix = branch_suffix,
                       verbose = verbose)
    }
  )

  #Add suffix to all target names within dynamic grouping variables
  if ('dynamic' %in% names(plan_modif)) {
    plan_modif$dynamic <- lapply(
      plan_modif$dynamic, function(command) {
        #print(command)
        modify_call_args(call=command,
                         in_regex = target_regex,
                         in_suffix = branch_suffix,
                         verbose = verbose)
      }
    )

  }


  return(plan_modif)
}
