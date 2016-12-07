\name{predict.DSA}
\alias{predict.DSA}
\title{Produces predicted values based on a DSA object}
\description{
  Predict returns predicted values based on a DSA object and potentially
  new data. If newdata is missing then predict produces the fitted
  values, otherwise the new data is used to produce predicted response
  values based on the model.selected from a call to the DSA routine.
}
\usage{
predict.DSA(object, newdata = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object returned from a call to the DSA routine.}
  \item{newdata}{An optional argument specifiying a new data frame to be
  used to obtain the fitted values.}
  \item{\dots}{ Currently ignored }
}
\value{
  The fitted values, or the results of calling glm.predict.
}
\keyword{models}% at least one, from doc/KEYWORDS

