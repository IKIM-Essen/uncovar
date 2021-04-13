{
    "pangolin strain (#SNPs)": function format(value) {
    $(function () {
        $('[data-toggle="popover"]').popover()
    })

    $(function () {
        $('[data-toggle="tooltip"]').tooltip()
    })

    if (value !== "no strain called") {
        var split = value.split(' ');
        var lineage = split[0];
        var link = `<a data-toggle="tooltip" data-placement="top" title="Linkout to cov-lineages" href='https://cov-lineages.org/lineages/lineage_${lineage}.html' target='_blank'>${lineage}</a>`;

        const table = `<div style="height: 200px; overflow-y: auto; white-space:pre-wrap;"><table class="table">
              <thead>
                <tr>
                  <th scope="col">Gene</th>
                  <th scope="col">Alteration</th>
                  <th scope="col">Present</th>
                </tr>
              </thead>
            <tbody>`;

        const table_end = "</tbody></table><div>";

        let inner_table = [];

        let cont = 0;
        let not_cont = 0;

        for (i = 1; i < split.length; i++) {
            let splitted_variant = split[i].split(":");
            let gene = splitted_variant[0];
            let alteration = splitted_variant[1];

            let contained = "";
            if (splitted_variant[2] === "true") {
                contained = "&#10003;" // Häkchen
                cont += 1;
            } else {
                contained = "&#10799;" // Kreuz
                not_cont += 1;
            }

            let row = `<tr><td scope="col">${gene}</td><td>${alteration}</td><td>${contained}</td></tr>`;
            inner_table.push(row);
        }

        let final_table = table + inner_table.join("") + table_end;
        let sum = cont + not_cont;

        let overview = `<a tabindex="0" role="button" href="#" data-toggle="popover" data-trigger="focus" data-html='true' title='Overview for ${lineage}' data-content='${final_table}'>(${cont}/${sum})</a>`

        if (split.length > 1) {
            return `${link} ${overview}`;
        } else {
            return `${link}`;
        }

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
                x = x + `<a tabindex="0" role="button" href="#" data-toggle="popover" data-trigger="focus" data-html='true' title='Gene: <a data-html="true" data-toggle="tooltip" data-placement="bottom" title="Linkout to gene in Ensembl genome browser" href="https://covid-19.ensembl.org/Sars_cov_2/Gene/Summary?g=${g}" target="_blank">${g}:${vats[g][i]}</a>' data-content='${tables[g]}'>${g}:${vats[g][i]}</a>`;
            }
        } else {
            var x = `<a tabindex="0" role="button" href="#" data-toggle="popover" data-trigger="focus" data-html='true' title='Gene: <a data-html="true" data-toggle="tooltip" data-placement="bottom" title="Linkout to gene in Ensembl genome browser" href="https://covid-19.ensembl.org/Sars_cov_2/Gene/Summary?g=${g}" target="_blank">${g}</a>' data-content='${tables[g]}'>${g}</a>`;
        }

        result.push(x);
    }
    
    result = result.join(", ");

    if (value.trim() !== "") {
        return result;
    } else {
        return "";
    }

},
    "variants of interest": function format(value) {
    return this["variant helper"](value, true);
},
    "other variants": function format(value) {
    let result = this["variant helper"](value, false);
    return result;
}
}
