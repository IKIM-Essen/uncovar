{
    "Pangolin Strain": function format(value) {
        $(function () {
            $('[data-toggle="tooltip"]').tooltip()
        })

        if (value !== "no strain called") {
            var link = `<a data-toggle="tooltip" data-placement="top" title="View ${value} on outbreak.info" href='https://outbreak.info/situation-reports?pango=${value}' target='_blank'>${value}</a>`;
            return link;
        } else {
            return value;
        }
    },
    "variant helper": function format(value, voi) {
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
            vafs[split[0]].push(split[2]);
            vats[split[0]].push(split[1]);
        }

        const table = `<table class="table">
                <thead>
                    <tr>
                    <th scope="col">Variant</th>
                    <th scope="col">VAF</th>
                    </tr>
                </thead>
                <tbody>`;

        const table_end = "</tbody></table>";

        var tables = {};

        for (g of genes) {
            var body = "";
            for (let i = 0; i < vafs[g].length; i++) {
                var row = `<tr><td scope="col">${vats[g][i]}</td><td>${vafs[g][i]}</td></tr>`;
                body = body + row;
            }
            tables[g] = table + body + table_end;
        }

        let result = [];
        let unique_genes = [...new Set(genes)];

        for (g of unique_genes) {
            if (voi) {
                var x = "";
                for (let i = 0; i < vats[g].length; i++) {
                    var x = `<a href="#" data-toggle="popover" data-trigger="focus" data-html='true' title='Gene: <a data-html="true" data-toggle="tooltip" data-placement="bottom" title="Linkout to gene in Ensembl genome browser" href="https://covid-19.ensembl.org/Sars_cov_2/Gene/Summary?g=${g}" target="_blank">${g}:${vats[g][i]}</a>' data-content='${tables[g]}'>${g}:${vats[g][i]}</a>`;
                    result.push(x);
                }
            } else {
                var x = `<a href="#" data-toggle="popover" data-trigger="focus" data-html='true' title='Gene: <a data-html="true" data-toggle="tooltip" data-placement="bottom" title="Linkout to gene in Ensembl genome browser" href="https://covid-19.ensembl.org/Sars_cov_2/Gene/Summary?g=${g}" target="_blank">${g}</a>' data-content='${tables[g]}'>${g}</a>`;
                result.push(x);
            }
        }
        
        result = result.join(", ");

        if (value.trim() !== "") {
            return result;
        } else {
            return "";
        }

    },
    "Variants of Interest": function format(value) {
        return this["variant helper"](value, true);
    },
    "Other Variants": function format(value) {
        let result = this["variant helper"](value, false);
        return result;
    }
}
