HTMLWidgets.widget({

  name: 'vennr',
  type: 'output',

  factory: function(el, width, height) {
    
    // Div object information
    var el = el;
    
    // Div object id
    var id = "#"+el.id
    
    // Venn diagram object
    var chart = venn.VennDiagram()
                    .width(width)
                    .height(height);

    return {

      renderValue: function(x) {
        
        // Set data in json format
        var sets = x.data;

        // Fill chart
        var div = d3.select(el);
        div.datum(sets).call(chart);

        // Area styling
        d3.selectAll(id+" .venn-circle path")
          .style("fill-opacity", 0.2)
          .style("stroke-width", 2)
          .style("stroke-opacity", 0.2)
          .style("stroke", "#444");

        // Text styling
        div.selectAll("text")
           .style("fill", "black")
           .style("font-size", x.settings.fontSize)
           .style("font-family", x.settings.fontFamily);
      },
      
      resize: function(width, height) {
        // console.log(height)
        // console.log(width)
        // Difficult to automatically resize venn.js charts
      }

    };
  }
});