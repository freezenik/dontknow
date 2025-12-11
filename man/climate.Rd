\name{climate}
\alias{climate}
\docType{data}

\title{Climate change attitudes with \dQuote{don't know} option from Eurobarometer 101.1}

\description{
  A processed subset of the Eurobarometer 101.1 (2024) survey containing
  an item on attitudes towards climate change together with selected
  socio-demographic covariates. The response is a Likert-type variable
  with an explicit \dQuote{don't know} category, used in the empirical
  application of the \code{\link{DK}} family.
}

\usage{
data("climate")
}

\format{
  A data frame with \emph{approximately} 480 observations (after
  filtering and listwise deletion of missing values) on the following
  9 variables:
  \describe{
    \item{climate}{
      Factor with five levels:
      \code{"don't know"},
      \code{"totally disagree"},
      \code{"tend to disagree"},
      \code{"tend to agree"},
      \code{"totally agree"}.
      This is based on Eurobarometer item \code{qb4\_7}:
      \emph{``Adapting to adverse climate change impacts can have positive
      outcomes for citizens in the EU''}. For model fitting with the
      \code{DK} family, this factor is transformed into a bivariate
      response \code{Y} as illustrated in the Examples.
    }
    \item{gender}{
      Factor indicating respondent's gender, with levels
      \code{"female"} and \code{"male"} (recoded from item \code{d10}).
    }
    \item{age}{
      Integer giving the respondent's age in years
      (recoded from item \code{d11}). Only respondents with
      \code{age >= 20} are retained.
    }
    \item{education}{
      Factor describing the highest level of education attained, with
      levels \code{"low"}, \code{"medium"}, and \code{"high"}.
      These categories are aggregated from the original education
      item \code{sd3b}:
      \code{"low"} combines pre-primary, primary, and lower secondary;
      \code{"medium"} combines upper secondary and post-secondary
      non-tertiary;
      \code{"high"} combines short-cycle tertiary, bachelor, master,
      and doctoral degrees.
      Observations with education \dQuote{completed abroad} or with
      \dQuote{don't know}/\dQuote{refusal} are excluded.
    }
    \item{occupation}{
      Factor indicating employment status with levels
      \code{"employed"}, \code{"nowork"}, and \code{"self"}
      (recoded from \code{d15a\_r1}). The reference level is
      \code{"self"}.
    }
    \item{hhs}{
      Factor giving the household size (number of persons aged
      at least one year) with levels
      \code{"one"}, \code{"two"}, \code{"three"}, and \code{"four+"}
      (recoded from \code{d40abc\_r}).
    }
    \item{urban}{
      Factor indicating degree of urbanization with levels
      \code{"rural"}, \code{"small"}, and \code{"large"},
      derived from item \code{d25}.
    }
    \item{country}{
      Factor with ISO country codes identifying the respondent's
      country of residence (from \code{isocntry}). East and West
      Germany are collapsed into a single category \code{"DE"}.
      The reference category is \code{"CY"} (Cyprus).
    }
    \item{politics}{
      Factor representing political interest with ordered levels
      \code{"not"}, \code{"low"}, \code{"medium"}, and \code{"strong"},
      recoded from item \code{polintr}. Responses \dQuote{don't know}
      or \dQuote{refusal} are treated as missing.
    }
  }
  All observations with missing values in any of the above variables are
  removed.
}

\details{
  The original Eurobarometer 101.1 micro data were imported using
  \pkg{haven} from the GESIS-provided Stata file
  \code{ZA7781\_v2-0-0.dta}. Categorical variables were reconstructed
  from the value labels, and a set of derived covariates was created to
  obtain a relatively simple, yet substantively meaningful, modelling
  setup.

  For use with the \code{\link{DK}} family, the univariate factor
  \code{climate} is transformed into a bivariate response
  \code{Y = (Y_1, Y_2)} where \code{Y_1} indicates whether the
  \dQuote{don't know} option is chosen and \code{Y_2} encodes the
  non-\dQuote{don't know} ordered categories. The exact transformation
  steps are given in the Examples and follow the implementation in the
  empirical application.
}

\source{
  European Commission and European Parliament, Brussels (2025).
  \emph{Eurobarometer 101.1 (2024)}.
  (ZA8841; Version 1.0.0) [Data set]. GESIS, Koeln.
  \doi{10.4232/1.14461}.
  \url{https://search.gesis.org/research_data/ZA7781}
}

\seealso{
  \code{\link{DK}}
}

\examples{\dontrun{data("climate")

## inspect response
table(climate$climate)

## construct the bivariate DK response
d <- climate

## create the response,
## reorder climate so DK is last.
d$climate2 <- factor(d$climate,
  levels = c("totally disagree", "tend to disagree", "tend to agree",
    "totally agree", "don't know" ))

d$y1 <- as.integer(d$climate == "don't know") ## 1 = DK, 0 = non-DK

## create y2 only for non-DK observations; arbitrary for DK (not used in LL).
tmp <- d$climate2
tmp[d$climate2 == "don't know"] <- NA
d$y2 <- as.integer(droplevels(tmp)) - 1L ## now 0, ..., 3 for 4 non-DK cats
d$y2[is.na(d$y2)] <- 0L                  ## any value for DK rows
d$Y <- cbind(d$y1, d$y2)

## fit a simple DK model
f <- Y ~ -1 + s(age) + la(gender, type = 2) + la(education, type = 3) |  ## y1
         -1 + la(occupation, type = 4) + la(country, type = 4) |         ## y2
         1                                                               ## rho

m <- gamlss2(f, family = DK(4), data = d)
summary(m)
}}

