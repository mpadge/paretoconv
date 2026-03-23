# Calculate Kummer's confluent hypergeometric function for real inputs
kummerM <- function(x, a, b) {
    .Call("chfm_call",
        as.double(x),
        as.double(a),
        as.double(b),
        PACKAGE = "paretoconv")
}
