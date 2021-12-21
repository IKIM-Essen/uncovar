{
    "Signatures": function(signature) {
        // taken from https://github.com/hodcroftlab/covariants/blob/338e94f7af9a7a1871434cfbb62af3f8d92ea90f/web/src/components/Common/parseAminoacidMutation.ts#L20
        if(signature == "Jaccard") {
            return signature
        }
        if(signature.includes(":")) {
            console.log(signature)
            const match = /^((?<gene>.*):)?(?<left>[*.a-z-]{0,1})(?<pos>(\d)*)(?<right>[*.a-z-]{0,1})$/i.exec(signature)
            console.log(match.groups)
            var gene = match.groups?.gene
            var ref = match.groups?.left
            var alt = match.groups?.right
            var pos = match.groups?.pos

            var gene_colors = {
                E: '#60AA9E',
                M: '#D9A456',
                N: '#D9776C',
                ORF10: '#E67030',
                ORF14: '#8EBC66',
                ORF1a: '#E59637',
                ORF1b: '#AABD52',
                ORF3a: '#C9957B',
                ORF6: '#5097BA',
                ORF7a: '#C4B945',
                ORF7b: '#75B681',
                ORF8: '#55AA81',
                ORF9b: '#D9AD3D',
                S: '#5097BA',
            }

            var aa_colors = {
                'A': '#EAEABA',
                'V': '#EAEA9F',
                'L': '#E1E177',
                'I': '#C9C94D',
                'B': '#AAAAAA',
                'C': '#E3F9B0',
                'D': '#E98F6D',
                'E': '#F7B080',
                'F': '#C7C88D',
                'G': '#C0C0C0',
                'H': '#D6F6FA',
                'K': '#CEC0F3',
                'M': '#C3ED3C',
                'N': '#F29290',
                'P': '#D2D1F8',
                'Q': '#F8C4E3',
                'R': '#A6ACEF',
                'S': '#D8B9D4',
                'T': '#F0D6E3',
                'W': '#86B0CC',
                'X': '#AAAAAA',
                'Y': '#8FC7D1',
                'Z': '#AAAAAA',
                '*': '#AAAAAA',
                '-': '#AAAAAA',
            }

            return `<span class="badge text-white p-0 m-1"><span class="p-1 rounded-left" style="background:${gene_colors[gene]};">${gene}:</span><span class="p-1" style="background:${aa_colors[ref]};">${ref}</span><span class="bg-secondary p-1">${pos}</span><span class="p-1 rounded-right" style="background:${aa_colors[alt]};">${alt}</span></span>`
        } else {
            return signature
        }
    },
    "VAF": function(vaf) {
        var lighting = (0.9 - parseFloat(vaf) * 0.4) * 100;
        return `<span class="badge p-0 m-1"><span class="p-1 rounded" style="background:hsl(190, 92%, ${lighting}%);">${vaf}</span></span>`
    },
    "Prob_present": function(prob) {
        var lighting = (0.9 - parseFloat(prob) * 0.4) * 100;
        return `<span class="badge p-0 m-1"><span class="p-1 rounded" style="background:hsl(138, 72%, ${lighting}%);">${prob}</span></span>`
    }
}
