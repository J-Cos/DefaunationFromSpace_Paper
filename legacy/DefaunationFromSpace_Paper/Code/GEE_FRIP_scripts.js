// =============================================================================
// LEGACY GEE FRIP SCRIPTS (JavaScript - Code Editor)
// =============================================================================
// These scripts were used in the original DefaunationFromSpace paper.
// Saved here for reference. Contains 4 scripts separated by comments.

// =============================================================================
// SCRIPT 1: Export base variable stack at MODIS native resolution
// =============================================================================

// load basin outlines
var hydro2 = ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_2');
var basins = hydro2
  .filter(
    ee.Filter.or(
      ee.Filter.eq('HYBAS_ID', 6020006540),
      ee.Filter.eq('HYBAS_ID', 1020018110)));
      
Map.addLayer(basins)

var getIEVI = function(year){
  var annualCollection = dataset
    .filter(ee.Filter.calendarRange(year, year, 'year'));
  var IEVI = annualCollection.mean();
  return(IEVI.set({'year': year}));
};

var yearList = ee.List.sequence(2001, 2023);

var dataset = ee.ImageCollection(ee.ImageCollection("MODIS/061/MOD17A3HGF").select("Npp").toList(23));

var modisProjection = ee.Image(dataset.first()).projection();

// IEVI
var annualIEVIList = yearList.map(getIEVI);
var medianIEVI = ee.ImageCollection.fromImages(annualIEVIList).median()

// flood frequency
var getReturnGroups = function(period){
  var periodCollection = ee.ImageCollection('JRC/CEMS_GLOFAS/FloodHazard/v1')
    .filterBounds(basins)
    .filter(ee.Filter.eq("return_period", period))
    .mosaic()
    .gte(0)
  return(periodCollection.set({'ReturnPeriod': period}));
};

var periodList = ee.List([10, 20, 50, 75, 100, 200, 500]);

var L = periodList.map(getReturnGroups);

var floodProjection = ee.Image(ee.ImageCollection('JRC/CEMS_GLOFAS/FloodHazard/v1').first()).projection();
var floodFreq = ee.ImageCollection.fromImages(L).sum().unmask().setDefaultProjection({crs: floodProjection, }).mask(ee.Image('MERIT/Hydro/v1_0_1').select('hnd').gt(0))

// Compute flooding per MODIS pixel.
var floodProb = floodFreq
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 65536
    });

var Forest = ee.ImageCollection('projects/JRC/TMF/v1_2023/AnnualChanges')
    .filterBounds(basins)
    .mosaic()
    .select('Dec2023')
    .eq(1)
    .setDefaultProjection(ee.Image(ee.ImageCollection('projects/JRC/TMF/v1_2023/AnnualChanges').first()).projection())
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 65536
    });

// median NPP
var variablesUnmasked= floodProb
  .addBands(medianIEVI)
  .addBands(Forest)
  .reproject({
    crs: modisProjection.crs(),
    crsTransform : modisProjection.getInfo()['transform']
  });

var vars = variablesUnmasked
  .updateMask(variablesUnmasked.select("Dec2023").gte(0.95))

//annual NPP
var variablesUnmasked_annual= floodProb
  .addBands(dataset
    .toBands()
    .rename(
      ["NPP01",
      "NPP02", "NPP03", "NPP04", "NPP05", "NPP06",
      "NPP07", "NPP08", "NPP09", "NPP10", "NPP11",
      "NPP12", "NPP13", "NPP14", "NPP15", "NPP16",
      "NPP17", "NPP18", "NPP19", "NPP20", "NPP21",
      "NPP22", "NPP23"
      ])
    .toDouble())
  .addBands(Forest)
  .reproject({
    crs: modisProjection.crs(),
    crsTransform : modisProjection.getInfo()['transform']
  });

var vars_annual = variablesUnmasked_annual
  .updateMask(variablesUnmasked_annual.select("Dec2023").gte(0.95))

// export images
Export.image.toAsset({
  image: vars,
  description: 'FloodingVars',
  assetId: 'FloodingVars',
  scale: 463.3127165279165,
  region: basins,
  maxPixels: 193675354
});

Export.image.toAsset({
  image: vars_annual,
  description: 'FloodingVarsAnnual',
  assetId: 'FloodingVarsAnnual',
  scale: 463.3127165279165,
  region: basins,
  maxPixels: 193675354
});


// =============================================================================
// SCRIPT 2: Export FRIP at multiple scales (from pre-computed base stack)
// =============================================================================

// load basin outlines
var hydro2 = ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_2');
var basins = hydro2
  .filter(
    ee.Filter.or(
      ee.Filter.eq('HYBAS_ID', 6020006540),
      ee.Filter.eq('HYBAS_ID', 1020018110)));

// function to export FRIP at multiple scales
var exportFRIPimages = function(scaleVal){
  
  var scale = ee.Number(scaleVal);
  
  var corr = vars. select("depth", "Npp")
      .reduceResolution({
        reducer: ee.Reducer.spearmansCorrelation(),
        maxPixels: 65536
      })
    .reproject({
      crs: vars.projection().crs(),
      scale: scale
    });
  
  var corrMask = corr.updateMask(corr.mask().gt(0.1))

  Export.image.toDrive({
    image: corrMask,
    description: 'FRIP_'+ scaleVal,
    folder: "GEE_FRIP",
    scale: scale,
    region: basins.geometry()
  });
}

// function to export annual FRIP at multiple scales
var exportAnnualFRIPimages = function(scaleVal){

  var scale = ee.Number(scaleVal);

  var NppYears =  ee.List(["NPP01",
        "NPP02", "NPP03", "NPP04", "NPP05", "NPP06",
        "NPP07", "NPP08", "NPP09", "NPP10", "NPP11",
        "NPP12", "NPP13", "NPP14", "NPP15", "NPP16",
        "NPP17", "NPP18", "NPP19", "NPP20", "NPP21",
        "NPP22", "NPP23"
        ])
  var yearList = ee.List.sequence(2001, 2023);
  
  var getAnnualCorrelation = function(index){
    var yearIm = varsannual
      .select([NppYears.get(index), "depth"])
    var corrAnnual = yearIm
        .reduceResolution({
          reducer: ee.Reducer.spearmansCorrelation(),
          maxPixels: 65536
        })
      .reproject({
        crs: varsannual.projection().crs(),
        scale: scale
      })
      .select("correlation")
      .set({'year': yearList.get(index)});
      
      var corrAnnualMask = corrAnnual.updateMask(corrAnnual.mask().gt(0.1))

    return(corrAnnualMask);
  }
  
  var CorrAnnualList = ee.List.sequence(0, 22).map(getAnnualCorrelation);
  var CorrAnnual = ee.ImageCollection.fromImages(CorrAnnualList)

  var corrIm = CorrAnnual.toBands()

  Export.image.toDrive({
    image: corrIm,
    description: 'FRIP_Annual_'+ scaleVal,
    folder: "GEE_FRIP",
    scale: scale,
    region: basins.geometry()
  });
}

// Define your list of scales
var scales = ee.List.sequence(5000, 100000, 5000);

scales.evaluate(function(scaleList) {
  scaleList.forEach(function(scale) {
    exportFRIPimages(scale);
  });
});

scales.evaluate(function(scaleList) {
  scaleList.forEach(function(scale) {
    exportAnnualFRIPimages(scale);
  });
});


// =============================================================================
// SCRIPT 3: Export masked predictor stack at MODIS native resolution
// =============================================================================

// load basin outlines
var hydro2 = ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_2');
var basins = hydro2
  .filter(
    ee.Filter.or(
      ee.Filter.eq('HYBAS_ID', 6020006540),
      ee.Filter.eq('HYBAS_ID', 1020018110)));

// Load MODIS MOD09A1
var modis = ee.ImageCollection('MODIS/061/MOD09A1')
  .select([
    'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03',
    'sur_refl_b04', 'sur_refl_b05', 'sur_refl_b06', 'sur_refl_b07'
  ]);

var getAnnuals = function(year){
  var annualCollection = modis
    .filter(ee.Filter.calendarRange(year, year, 'year'));
  var annuals = annualCollection.mean();
  return(annuals.set({'year': year}));
};

var yearList = ee.List.sequence(2001, 2023);

var annualList = yearList.map(getAnnuals);
var medianModis = ee.ImageCollection.fromImages(annualList).median()

// NDVI
var ndvi = medianModis.normalizedDifference(['sur_refl_b02', 'sur_refl_b01']).rename('ndvi');
medianModis = medianModis.addBands(ndvi);

// flood frequency
var getReturnGroups = function(period){
  var periodCollection = ee.ImageCollection('JRC/CEMS_GLOFAS/FloodHazard/v1')
    .filterBounds(basins)
    .filter(ee.Filter.eq("return_period", period))
    .mosaic()
    .gte(0)
  return(periodCollection.set({'ReturnPeriod': period}));
};

var periodList = ee.List([10, 20, 50, 75, 100, 200, 500]);

var L = periodList.map(getReturnGroups);

var floodProjection = ee.Image(ee.ImageCollection('JRC/CEMS_GLOFAS/FloodHazard/v1').first()).projection();
var floodFreq = ee.ImageCollection.fromImages(L).sum().unmask().setDefaultProjection({crs: floodProjection, }).mask(ee.Image('MERIT/Hydro/v1_0_1').select('hnd').gt(0))

var floodProb = floodFreq
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 65536
    })
    .reproject({
    crs: modis.first().projection().crs(),
    crsTransform : modis.first().projection().getInfo()['transform']
    });

var elev = ee.Image('CGIAR/SRTM90_V4').select('elevation')
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 65536
    })
    .reproject({
    crs: modis.first().projection().crs(),
    crsTransform : modis.first().projection().getInfo()['transform']
    });

var HAND = ee.Image('MERIT/Hydro/v1_0_1').select('hnd')
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 65536
    })
    .reproject({
    crs: modis.first().projection().crs(),
    crsTransform : modis.first().projection().getInfo()['transform']
    });

var Forest = ee.ImageCollection('projects/JRC/TMF/v1_2023/AnnualChanges')
    .filterBounds(basins)
    .mosaic()
    .select('Dec2023')
    .eq(1)
    .setDefaultProjection(ee.Image(ee.ImageCollection('projects/JRC/TMF/v1_2023/AnnualChanges').first()).projection())
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 65536
    })
    .reproject({
      crs: modis.first().projection().crs(),
      crsTransform : modis.first().projection().getInfo()['transform']
    });

var maskedPredictors = medianModis
  .setDefaultProjection(modis.first().projection())
  .addBands(floodProb)
  .addBands(elev)
  .addBands(HAND)
  .updateMask(Forest.select("Dec2023").gte(0.95))

Export.image.toAsset({
  image: maskedPredictors,
  description: 'maskedPredictors',
  assetId: 'maskedPredictors',
  scale: 463.3127165279165,
  region: basins,
  maxPixels: 1e13
});


// =============================================================================
// SCRIPT 4: Export predictor stack at multiple scales
// =============================================================================

// load basin outlines
var hydro2 = ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_2');
var basins = hydro2
  .filter(
    ee.Filter.or(
      ee.Filter.eq('HYBAS_ID', 6020006540),
      ee.Filter.eq('HYBAS_ID', 1020018110)));

var exportMODISimages = function(scaleVal){
  
  var scale = ee.Number(scaleVal);

  var ModisStack = im
    .reduceResolution({
      reducer: ee.Reducer.mean().combine({
        reducer2: ee.Reducer.variance(),
        sharedInputs: true
      }),
      maxPixels: 65536
    })
    .reproject({
      crs: im.projection().crs(),
      scale: scale
    })

  var ModisStack = ModisStack.addBands(ModisStack.select("ndvi_mean").mask().rename("LocalCover"))
    .toFloat();

  Export.image.toDrive({
    image: ModisStack,
    description:  'PredictorStack_'+ scaleVal,
    scale: scale,
    region: basins.geometry(),
    folder: "GEE_PredictorStack",
    maxPixels: 1e13
  });
}

var scales = ee.List.sequence(5000, 100000, 5000);

scales.evaluate(function(scaleList) {
  scaleList.forEach(function(scale) {
    exportMODISimages(scale);
  });
});
