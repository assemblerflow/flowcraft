class Params {

    static void check(Map params) {

        // Checks genomeSize for type
        try {
            params.genomeSize as Double
        } catch (e) {
            print_error("The genomeSize option must be a number")
        }

        // Checks minCoverage for type
        try {
            params.minCoverage as Double
        } catch (e) {
            print_error("the minCoverage option must be a number")
        }

        // Check if fastqc adapters file exists
        if (!params.adapters.equalsIgnoreCase("none")) {
            File f = new File(params.adapters)
            if (!f.exists()) {
                print_error("The provided adapters file does " +
                            "not exist ($params.adapters)")
            }
        }

        // Check for trimmomatic parameters
        try {
            params.trimLeading as Double
            params.trimTrailing as Double
            params.trimMinLength as Double
        } catch (e) {
            print_error("The trimLeading ($params.trimLeading), " +
                        "trimTrailing ($params.trimTrailing) and " +
                        "trimMinLength ($params.trimMinLength) " +
                        "options must be numbers")
        }

        // Check for Spades parameters
        [
            "spadesMincoverage": params.spadesMinCoverage,
            "spadesMinKmerCoverage": params.spadesMinKmerCoverage,
            "spadesMinContigLen": params.spadesMinContigLen,
            "spadesMaxContigs": params.spadesMaxContigs
        ].each { k, v ->
            try {
                v as Integer
            } catch (e) {
                print_error("The spades parameter $k ($v) must be an integer")
            }
         }

    }

    static def print_error(String msg) {

        println "\nERROR: $msg"
        System.exit(1)

    }

}