{
    "pangolin strain (#SNPs)": function format(value) {
        $(function () {
            $('[data-toggle="tooltip"]').tooltip()
        })

        if (value !== "no strain called") {
            var lineage = value.split(' ')[0];
            var link = "<a data-toggle=\"tooltip\" data-placement=\"top\" title=\"Linkout to cov-lineages\" href='https://cov-lineages.org/lineages/lineage_" + lineage + ".html' target='_blank'>" + lineage + "</a>";
            return link + value.split(lineage).pop();
        } else {
            return value;3
        }
    },
    "variants of interest": function format(value) {
        $(function () {
            $('[data-toggle="popover"]').popover()
        })

        $(function () {
            $('[data-toggle="tooltip"]').tooltip()
        })

        $('.popover-dismiss').popover({
            trigger: 'focus'
        })

        var variants = value.split(' ');
        var genes = [];
        var vafs = {};
        var vats = {};
        for (var v of variants) {
            var split = v.split(':');
            genes.push(split[0]);
        }

        for (var g of genes) {
            vafs[g] = [];
            vats[g] = [];
        }

        for (var v of variants) {
            var split = v.split(':');
            vafs[split[0]].push(split[1]);
            vats[split[0]].push(split[2]);
        }

        const table = "<table class=\"table\">\n" +
            "  <thead>\n" +
            "    <tr>\n" +
            "      <th scope=\"col\">Variant</th>\n" +
            "      <th scope=\"col\">VAF</th>\n" +
            "    </tr>\n" +
            "  </thead>" +
            "<tbody>";

        const table_end = "</tbody></table>";

        var tables = {};

        for (g of genes) {
            var body = "";
            for (let i = 0; i < vafs[g].length; i++) {
                var row = "    <tr><td scope=\"col\">" + vafs[g][i] + "</td><td>"  + vats[g][i] + "</td></tr>";
                body = body + row;
            }
            tables[g] = table + body + table_end;
        }

        let result = "";
        let unique_genes = [...new Set(genes)];

        for (g of unique_genes) {
            var x = "<a tabindex=\"0\" class=\"btn btn-link\" data-toggle=\"popover\" data-trigger=\"focus\" data-html='true' title='Gene: <a data-html=\"true\" data-toggle=\"tooltip\" data-placement=\"bottom\" title=\"Linkout to gene in Ensembl genome browser\" href=\"https://covid-19.ensembl.org/Sars_cov_2/Gene/Summary?g=" + g + "\" target=\"_blank\">" + g + "</a>' data-content='" + tables[g] + "'>" + g +"</a>";
            result = result + x + " ";
        }

        if (value.trim() !== "") {
            return result;
        } else {
            return "";
        }

        }
}
