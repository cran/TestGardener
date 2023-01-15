\name{NEWS}
\title{Recent changes to the TestGardner package}
\section{Changes in TestGardener version 3.1.3 (2022-20-12)}{
\itemize{
  \item{We observe that the function used to define the best estimate of 
    the score index for an examinee about 10 to 15 percent of the time
    has large nearly flat intervals of low values.  Moreover, it is
    often the case that there are multiple minimizing locations, and
    function thetafun does not always identify the lowest of these.
    New function thetasearch uses a fast grid search to identify the
    lowest point, and replaces the location identified by thetafun where
    required by the integer closest to the global minimum.  This is 
    usually done after the cycles are completed as a final step.
    }
  }
}
\section{Changes in TestGardener version 3.1.2 (2022-20-12)}{
\itemize{
  \item{The test analyst no longer needs to construct a "grbg"" object in order
    handle missing of illegitimate responses by treating them as an
    additional option, the garbage option as we call it.  This is now
    done automatically within function make.dataList.  Function make.dataList
    no longer needs an argument "grbg"."
    In order to makethis work, the user does need to adhere to a discipline.  
    For each item, the labels or indices of its options must be the first 
    positive integers.  
    For the quantitative SweSAT test, for example, 68 of the items have option 
    indices 1,...,4, and 12 have indices 1,...,5. For this large test given
    to over 55,000 examinees, the garbage option is essential for all items.
    It is assigned within make.dataList to an option with the index one larger 
    than required for the legitimate options, so that after make_dataList
    has done its work, the 68 items have indices 1,...,5 and the others
    indices 1,...,6.  
    Any integer within the data matrix U that exceeds the highest legitimate 
    index is re-assigned to the corresponding garbage index.  
    Do be aware that this process can change the U matrix that you set up 
    initially.
    Users often want to have non-integer labels for options.  This is still
    available by setting these up in the "optList" list vector, and these
    will appear in the most of the plots.
    }
  \item{This version of TestGardener also has various code improvements
    that speed up computation and improve the information in plots.
    }
  }
}
\section{Changes in TestGardener version 3.0.1 (2022-22-07)}{
\itemize{
  \item{Plots using the score index as the abscissa now use arc length along
    along the test information curve as the abscissa, since this interval is
    invariant under warpings of the index and therefore on a ratio scale. The 
    codes in the demo and vignette folders have been modified accordingly.
    }
  }
}
