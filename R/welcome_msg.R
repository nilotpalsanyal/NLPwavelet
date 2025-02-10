.onAttach <- function(libname, pkgname) {
  message <- c("\n Welcome to NLPwavelet!",
               "\n \n Website: https://nilotpalsanyal.github.io/NLPwavelet/",
               "\n Bug report: https://github.com/nilotpalsanyal/NLPwavelet/issues")
  packageStartupMessage(message)
}