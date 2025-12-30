#' Normalize Chemical Name Capitalization
#'
#' Converts chemical names to proper title case while preserving stereochemistry
#' prefixes, abbreviations, and other chemically meaningful capitalization.
#'
#' @param name Character vector of chemical names to normalize
#'
#' @return Character vector with normalized capitalization
#' @export
#'
#' @examples
#' normalize_chem_name("3-HYDROXYBUTANOIC ACID")
#' # Returns: "3-Hydroxybutanoic acid"
#' normalize_chem_name("D-ASPARTATE")
#' # Returns: "D-Aspartate"
normalize_chem_name <- function(name) {
  
  # Return NA/empty as-is
  if (is.na(name) || name == "") return(name)
  
  # Quick check if name is already reasonably formatted
  # Look for names that are already in proper title case
  total_chars <- nchar(gsub("[^A-Za-z]", "", name))
  upper_chars <- nchar(gsub("[^A-Z]", "", name))
  lower_chars <- nchar(gsub("[^a-z]", "", name))
  
  if (total_chars > 0) {
    # Check if it's already in good title case format:
    # - Starts with capital letter
    # - Has lowercase letters 
    # - Internal words after spaces/hyphens are capitalized appropriately
    # - Not mostly uppercase (avoid all-caps names)
    upper_ratio <- upper_chars / total_chars
    
    # Must start with capital and not be mostly uppercase
    if (upper_ratio < 0.8 && grepl("^[A-Z]", name)) {
      # Split into words to check title case pattern
      words <- unlist(strsplit(name, "[ -]+"))
      words <- words[nchar(words) > 0]  # Remove empty
      
      # Check if most words follow title case (start with capital, rest lowercase)
      title_case_words <- sum(grepl("^[A-Z][a-z]*$|^[0-9]", words))
      total_words <- length(words)
      
      # If ALL words are in proper title case, consider it well-formatted
      if (total_words > 0 && title_case_words == total_words) {
        # Only apply critical fixes for parenthetical stereochemistry
        result <- gsub("\\(r\\)-", "(R)-", name)
        result <- gsub("\\(s\\)-", "(S)-", result)
        result <- gsub("\\(r,s\\)-", "(R,S)-", result)
        result <- gsub("\\(rs\\)-", "(RS)-", result)
        result <- gsub("\\(sr\\)-", "(SR)-", result)
        result <- gsub("\\(([0-9]+)h\\)-", "(\\1H)-", result)
        result <- gsub("\\(([+-])\\)-", "(\\1)-", result)
        return(result)
      }
    }
  }
  
  # Convert to lowercase first
  result <- tolower(name)
  
  # Title case: capitalize first letter of each word
  # Split by spaces and hyphens (but keep hyphens)
  words <- strsplit(result, "(?<=[ ])|(?=[ ])", perl = TRUE)[[1]]
  
  # Capitalize first letter of first word and words after spaces
  capitalize_word <- function(word) {
    if (nchar(word) == 0 || word == " ") return(word)
    # Don't capitalize if it's just a number or starts with number
    if (grepl("^[0-9]", word)) {
      # Find first letter position and capitalize it
      pos <- regexpr("[a-z]", word)
      if (pos > 0) {
        substr(word, pos, pos) <- toupper(substr(word, pos, pos))
      }
      return(word)
    }
    paste0(toupper(substr(word, 1, 1)), substr(word, 2, nchar(word)))
  }
  
  # Apply to first word and words after spaces
  for (i in seq_along(words)) {
    if (i == 1 || (i > 1 && words[i-1] == " ")) {
      words[i] <- capitalize_word(words[i])
    }
  }
  result <- paste0(words, collapse = "")
  
  # Stereochemistry prefixes should be specific case:
  # D-, L-, R-, S- should be UPPERCASE
  result <- gsub("\\b([dlrs])-", "\\U\\1-", result, perl = TRUE)
  
  # Greek letters and positional prefixes should be lowercase
  result <- gsub("\\bAlpha-", "alpha-", result)
  result <- gsub("\\bBeta-", "beta-", result)
  result <- gsub("\\bGamma-", "gamma-", result)
  result <- gsub("\\bDelta-", "delta-", result)
  result <- gsub("\\bEpsilon-", "epsilon-", result)
  result <- gsub("\\bOmega-", "omega-", result)
  result <- gsub("\\bMyo-", "myo-", result)
  result <- gsub("\\bChiro-", "chiro-", result)
  result <- gsub("\\bCis-", "cis-", result)
  result <- gsub("\\bTrans-", "trans-", result)
  result <- gsub("\\bSec-", "sec-", result)
  result <- gsub("\\bTert-", "tert-", result)
  result <- gsub("\\bN-", "N-", result)  # Keep N- uppercase (nitrogen)
  result <- gsub("\\bO-", "O-", result)  # Keep O- uppercase (oxygen)
  result <- gsub("\\bS-", "S-", result)  # Keep S- uppercase (sulfur) - but watch for S-stereochem
  result <- gsub("\\bM-", "m-", result)  # meta position lowercase
  result <- gsub("\\bP-", "p-", result)  # para position lowercase
  
  # Handle special cases
  # DL- should be uppercase
  result <- gsub("\\bDl-", "DL-", result)
  
  # Rac- for racemic
  result <- gsub("\\bRac-", "rac-", result)
  
  # Fix parenthetical stereochemistry and chemical notations (should be uppercase)
  result <- gsub("\\(r\\)-", "(R)-", result)
  result <- gsub("\\(s\\)-", "(S)-", result)
  result <- gsub("\\(r,s\\)-", "(R,S)-", result)
  result <- gsub("\\(rs\\)-", "(RS)-", result)
  result <- gsub("\\(sr\\)-", "(SR)-", result)
  result <- gsub("\\(([0-9]+)h\\)-", "(\\1H)-", result)  # (2h) -> (2H)
  result <- gsub("\\(([+-])\\)-", "(\\1)-", result)      # Handle (+)/(-) signs
  
  # Common chemical suffixes should be lowercase (when not first word)
  result <- gsub(" Acid\\b", " acid", result)
  result <- gsub(" Ester\\b", " ester", result)
  result <- gsub(" Aldehyde\\b", " aldehyde", result)
  result <- gsub(" Amine\\b", " amine", result)
  result <- gsub(" Alcohol\\b", " alcohol", result)
  result <- gsub(" Oxide\\b", " oxide", result)
  result <- gsub(" Phosphate\\b", " phosphate", result)
  result <- gsub(" Sulfate\\b", " sulfate", result)
  result <- gsub(" Fatty\\b", " fatty", result)
  result <- gsub(" Salt\\b", " salt", result)
  result <- gsub(" Conjugate\\b", " conjugate", result)
  result <- gsub(" Peptide\\b", " peptide", result)
  result <- gsub(" Glycol\\b", " glycol", result)
  result <- gsub(" Succinate\\b", " succinate", result)
  result <- gsub(" Vanillate\\b", " vanillate", result)
  result <- gsub(" Monoacetate\\b", " monoacetate", result)
  result <- gsub(" Dimer\\b", " dimer", result)
  
  # Handle (R)- and (S)- stereochemistry in parentheses

  result <- gsub("\\(r\\)-", "(R)-", result)
  result <- gsub("\\(s\\)-", "(S)-", result)
  result <- gsub("\\(r,s\\)-", "(R,S)-", result)
  result <- gsub("\\(s,r\\)-", "(S,R)-", result)
  
  # Handle N(pai) type patterns - restore PAI to uppercase
  result <- gsub("N\\(pai\\)", "N(pai)", result, ignore.case = TRUE)
  
 # Fix common abbreviations (case-insensitive match, restore proper case)
  result <- gsub("\\bdihome\\b", "DiHOME", result, ignore.case = TRUE)
  result <- gsub("\\bhode\\b", "HODE", result, ignore.case = TRUE)
  result <- gsub("\\bhete\\b", "HETE", result, ignore.case = TRUE)
  result <- gsub("\\bcehc\\b", "CEHC", result, ignore.case = TRUE)
  result <- gsub("\\btmcp\\b", "TMCP", result, ignore.case = TRUE)
  result <- gsub("\\bcoa\\b", "CoA", result, ignore.case = TRUE)
  result <- gsub("\\bnadph\\b", "NADPH", result, ignore.case = TRUE)
  result <- gsub("\\bnadp\\b", "NADP", result, ignore.case = TRUE)
  result <- gsub("\\bnadh\\b", "NADH", result, ignore.case = TRUE)
  result <- gsub("\\bnad\\b", "NAD", result, ignore.case = TRUE)
  result <- gsub("\\batp\\b", "ATP", result, ignore.case = TRUE)
  result <- gsub("\\badp\\b", "ADP", result, ignore.case = TRUE)
  result <- gsub("\\bamp\\b", "AMP", result, ignore.case = TRUE)
  result <- gsub("\\bgtp\\b", "GTP", result, ignore.case = TRUE)
  result <- gsub("\\bgdp\\b", "GDP", result, ignore.case = TRUE)
  result <- gsub("\\bgmp\\b", "GMP", result, ignore.case = TRUE)
  result <- gsub("\\bdna\\b", "DNA", result, ignore.case = TRUE)
  result <- gsub("\\brna\\b", "RNA", result, ignore.case = TRUE)
  
  return(result)
}

#' Vectorized version of normalize_chem_name
#' @param names Character vector
#' @return Character vector with normalized names
#' @export
normalize_chem_names <- function(names) {
  sapply(names, normalize_chem_name, USE.NAMES = FALSE)
}
