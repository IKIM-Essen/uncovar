{
  Mutations: function (signature) {
    // prettier-ignore
    // taken from https://github.com/hodcroftlab/covariants/blob/338e94f7af9a7a1871434cfbb62af3f8d92ea90f/web/src/components/Common/parseAminoacidMutation.ts#L20
    if(signature == "Lineage") {
            return `<span class="badge p-0 m-1"><span class="p-1 rounded bg-secondary text-white">${"Lineage:"}</span></span>`
        }
    if (signature == "Similarity") {
      return `<span class="badge p-0 m-1"><span class="p-1 rounded bg-secondary text-white">${"Similarity:"}</span></span>`;
    }
    if (signature.includes(":")) {
      console.log(signature);
      const match =
        /^((?<gene>.*):)?(?<left>[*.a-z-]{0,1})(?<pos>(\d)*)(?<right>[*.a-z-]{0,1})$/i.exec(
          signature
        );
      console.log(match.groups);
      var gene = match.groups?.gene;
      var ref = match.groups?.left;
      var alt = match.groups?.right;
      var pos = match.groups?.pos;

      var gene_colors = {
        E: "#60AA9E",
        M: "#D9A456",
        N: "#D9776C",
        ORF10: "#E67030",
        ORF14: "#8EBC66",
        ORF1a: "#E59637",
        ORF1b: "#AABD52",
        ORF3a: "#C9957B",
        ORF6: "#5097BA",
        ORF7a: "#C4B945",
        ORF7b: "#75B681",
        ORF8: "#55AA81",
        ORF9b: "#D9AD3D",
        S: "#5097BA",
      };

      var aa_colors = {
        A: "#EAEABA",
        V: "#EAEA9F",
        L: "#E1E177",
        I: "#C9C94D",
        B: "#AAAAAA",
        C: "#E3F9B0",
        D: "#E98F6D",
        E: "#F7B080",
        F: "#C7C88D",
        G: "#C0C0C0",
        H: "#D6F6FA",
        K: "#CEC0F3",
        M: "#C3ED3C",
        N: "#F29290",
        P: "#D2D1F8",
        Q: "#F8C4E3",
        R: "#A6ACEF",
        S: "#D8B9D4",
        T: "#F0D6E3",
        W: "#86B0CC",
        X: "#AAAAAA",
        Y: "#8FC7D1",
        Z: "#AAAAAA",
        "*": "#AAAAAA",
        "-": "#AAAAAA",
      };

      return `<span class="badge text-white p-0 m-1"><span class="p-1 rounded-left" style="background:${gene_colors[gene]};">${gene}:</span><span class="p-1" style="background:${aa_colors[ref]};">${ref}</span><span class="bg-secondary p-1">${pos}</span><span class="p-1 rounded-right" style="background:${aa_colors[alt]};">${alt}</span></span>`;
    } else {
      var nucl_colors = {
        A: "#B6EE92",
        T: "#B4D3FA",
        G: "#FFD63F",
        C: "#FBBFAB",
        // 'U': '#8A89FF',
        // 'R': '#FFFE80',
        // 'Y': '#E180FF',
        // 'S': '#FF9B80',
        // 'W': '#80FFF2',
        // 'M': '#CE8834',
        // 'K': '#90B82C',
        // 'D': '#C7FFB9',
        // 'B': '#F8C1C0',
        // 'V': '#FFE3B9',
        // 'H': '#BFD8F9',
        N: "#AAAAAA",
        X: "#AAAAAA",
        "-": "#AAAAAA",
      };

      const match =
        /^(?<left>[.a-z-]{0,1})(?<pos>(\d)*)(?<right>[.a-z-]{0,1})$/i.exec(
          signature
        );

      var ref = match.groups?.left;
      var alt = match.groups?.right;
      var pos = match.groups?.pos;

      return `<span class="badge text-white p-0 m-1""><span class="p-1 rounded-left" style="background:${nucl_colors[ref]};">${ref}</span><span class="bg-secondary p-1";>${pos}</span><span class="p-1 rounded-right" style="background:${nucl_colors[alt]};">${alt}</span></span>`;
    }
  },
  Frequency: function (vaf) {
    var lighting = (0.9 - parseFloat(vaf) * 0.4) * 100;
    return `<span class="badge p-0 m-1"><span class="p-1 rounded" style="background:hsl(190, 92%, ${lighting}%);">${vaf}</span></span>`;
  },
  Probability: function (prob) {
    var lighting = (0.9 - parseFloat(prob) * 0.4) * 100;
    return `<span class="badge p-0 m-1"><span class="p-1 rounded" style="background:hsl(138, 72%, ${lighting}%);">${prob}</span></span>`;
  },
  "lineage helper": function (value) {
    if (isNaN(value)) {
      var variant_colors = {
        '19A': { bg: '#0F2B9D', fg: '#FFFFFF' },
        '19B': { bg: '#302278', fg: '#FFFFFF' },
        '20A': { bg: '#4A1A5A', fg: '#FFFFFF' },
        '20B': { bg: '#691136', fg: '#FFFFFF' },
        '20C': { bg: '#840917', fg: '#FFFFFF' },
        '20D': { bg: '#8E0837', fg: '#FFFFFF' },
        '20E': { bg: '#930876', fg: '#FFFFFF' },
        '20F': { bg: '#9609A7', fg: '#FFFFFF' },
        '20G': { bg: '#9B0AE6', fg: '#FFFFFF' },
        '20H': { bg: '#4D21AD', fg: '#FFFFFF' },
        '20I': { bg: '#403FC6', fg: '#FFFFFF' },
        '20J': { bg: '#3F63CF', fg: '#FFFFFF' },
        '21A': { bg: '#4783C8', fg: '#000000' },
        '21I': { bg: '#539BB5', fg: '#000000' },
        '21J': { bg: '#63AC9B', fg: '#000000' },
        '21B': { bg: '#77B67F', fg: '#000000' },
        '21C': { bg: '#8EBC66', fg: '#000000' },
        '21D': { bg: '#A8BD53', fg: '#000000' },
        '21E': { bg: '#C1BA47', fg: '#000000' },
        '21F': { bg: '#D6B03F', fg: '#000000' },
        '21G': { bg: '#E39C39', fg: '#000000' },
        '21H': { bg: '#E67C33', fg: '#000000' },
        '21K': { bg: '#E1532B', fg: '#000000' },
        '21L': { bg: '#DB2823', fg: '#000000' },
      };

      if (value == "x") {
        return `<span class="badge p-0 m-1"><span class="p-1 rounded" style="background:#24DC5B;">${"\u2713"}</span></span>`;
      } else {
        const match = /^(?<parent>.{3})\s(?<version>.+)\s?.*$/i.exec(value);

        var parent = match.groups?.parent;
        var version = match.groups?.version;

        return `<span class="badge p-0 m-1 text-white"><span class="p-1 rounded-left" style="background:${variant_colors[parent].bg}; color:${variant_colors[parent].fg}">${parent}</span><span class="bg-secondary p-1 rounded-right">${version}</span></span>`;
        // return this["variant helper"](value);
      }
    } else {
      if (value == "") {
        return value;
      } else {
        var lighting = (0.9 - parseFloat(value) * 0.4) * 100;
        value = parseFloat(value) * 100;
        value = value.toFixed(2)
        return `<span class="badge p-0 m-1"><span class="p-1 rounded-left" style="background:hsl(138, 72%, ${lighting}%);">${value}</span><span class="p-1 rounded-right" style="background:hsl(138, 72%, ${lighting}%);">${"%"}</span></span>`;
      }
    }
  },
  "Highest similarity": function (value) {
    let result = this["lineage helper"](value);
    return result;
  },
  "2nd": function (value) {
    let result = this["lineage helper"](value);
    return result;
  },
  "3rd": function (value) {
    let result = this["lineage helper"](value);
    return result;
  },
  "4th": function (value) {
    let result = this["lineage helper"](value);
    return result;
  },
  "5th": function (value) {
    let result = this["lineage helper"](value);
    return result;
  },
};
