/*
A KBase module: nkk_compHelloWorld
*/

module nkk_compHelloWorld {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_nkk_compHelloWorld(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
