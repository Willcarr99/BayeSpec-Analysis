## Help file

cat("\n###########################################\n")
cat("This is the help file\n\n")

cat("select.single()\n",
    "Select a single beak with no background\n",
    "Returns (channel, counts, energy)\n\n")

cat("select.triple()\n",
    "Select a peak with background sidebands\n",
    "Returns (channel, counts, background counts,\n",
    "         background fit uncertainties, energy)\n\n")
