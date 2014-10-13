 \name{newQSarray}
 \alias{newQSarray}
 \title{The qusage Array Object}
 \description{
   The constructor for the QSarray object. Should primarily be used internally by \code{\link{qusage}} or \code{\link{makeComparison}}. See \code{\link{QSarray-class}} for a full description of the fields in the QSarray object.
 }
 \usage{
newQSarray(obj=NULL,
        \dots
       )
 }
 \arguments{ 
    \item{obj}{The features of QSarray can be supplied as a list of objects. The objects in the list must be named appropriately. See \code{\link{QSarray-class}} for a description of the parameters which can be stored in the QSarray object.}
    \item{\dots}{The fields of the QSarray object can also be specified individually. If \var{obj} is specified and additional fields are provided, the parameters will be combined together into a single QSarray object, with the parameters specified by \dots replacing those in \var{obj} (this will also produce a warning).}

}
\details{
This is the constructor for use in creating QSarray objects. This is primarily intended for internal use, but advanced users may find it useful to construct their own QSarray objects without going through the process of running qusage.

In order to create a QSarray object from scratch, the constructor requires the following three fields: \code{mean}, \code{sd}, and \code{dof}. All other fields can be either left blank or added after. Note that in some cases, various methods will not be able to run without more information. For a complete list of the fields that the QSarray object can contain, refer to \code{\link{QSarray-class}}.
}





