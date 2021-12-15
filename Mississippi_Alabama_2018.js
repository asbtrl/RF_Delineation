//var points = SE_TF.merge(SE_Sea).merge(SE_BarLand).merge(SE_VegLand).merge(SE_UrbLand);
//var state1 = ee.FeatureCollection('TIGER/2016/States')
//    .filter(ee.Filter(
//        ee.Filter.eq('NAME', 'Louisiana')));
//var aoi1 = coast_area.intersection(state1, ee.ErrorMargin(1));
        
var state2 = ee.FeatureCollection('TIGER/2016/States')
    .filter(ee.Filter(
        ee.Filter.eq('NAME', 'Mississippi')));
var aoi2 = coast_area.intersection(state2, ee.ErrorMargin(1));

var state3 = ee.FeatureCollection('TIGER/2016/States')
    .filter(ee.Filter(
        ee.Filter.eq('NAME', 'Alabama')));
var aoi3 = coast_area.intersection(state3, ee.ErrorMargin(1));

//var region = aoi1.union(aoi2, ee.ErrorMargin(1)).union(aoi3, ee.ErrorMargin(1)).union(add_up, ee.ErrorMargin(1));
//var region = aoi1;
var region = aoi2.union(aoi3, ee.ErrorMargin(1));

// 0. Global Variables
var globOptions = {
  versionID: '_SR',
  outFolder: 'SR',
    startDate: '2018-01-01',
  endDate: '2018-12-31',
  bandSelect: ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'],
  bands8: ['B2', 'B3', 'B4', 'B5', 'B6', 'B7'],
  bands7: ['B1', 'B2', 'B3', 'B4', 'B5', 'B7'], 
  maskAltitude: 100,  
  maskDepth: -100, 
  maskDistance: 15000,
  maskApplySRTM: false,
  parallelScale: 8,
  trainingValidationRatio: 0.8,
  nTrees: 10, 
  outScale: 30, 
  conPixels: 100
};

var distFilter = ee.Filter.withinDistance({
  distance: 100000,
  leftField: '.geo',
  rightField: '.geo',
  maxError: 10
});

// Define a saveAll join.
var distSaveAll = ee.Join.saveAll({
  matchesKey: 'points',
  measureKey: 'distance'
});

// Apply the join.
var spatialJoined_points = distSaveAll.apply(points,region, distFilter);
print(spatialJoined_points);

var Regional_coast = US_coast.filterBounds(region);
var aoi=region;

// 1. Functions
var landsatFunctions = {

  applyNDWI: function(image) {
    // apply NDWI to an image
    var ndwi = image.normalizedDifference(['green','nir']);
    return ndwi.select([0], ['ndwi']);
  },

  applyMNDWI: function(image) {
    // apply MNDWI to an image
    var mndwi = image.normalizedDifference(['green','swir1']);
    return image.select([0], ['mndwi']);
  },

  applyAWEI: function(image) {
    // apply AWEI to an image
    var awei = image.expression("b('blue')+2.5*b('green')-1.5*(b('nir')+b('swir1'))-0.25*b('swir2')");
    return awei.select([0], ['awei']);
  },

  applyNDVI: function(image) {
    // apply NDVI to an image
    var ndvi = image.normalizedDifference(['nir','red']);
    return ndvi.select([0], ['ndvi']);
  },
  
  applyEVI: function(image) {
    // apply EVI to an image
    var evi = image.expression("2.5*(b('nir')-b('red'))/(b('nir')+6*b('red')-7.5*b('blue')+1)");
    return evi.select([0], ['evi']);
  },
  
  applyLSWI: function(image) {
    // apply LSWI to an image
    var lswi = image.normalizedDifference(['nir','swir1']);
    return lswi.select([0], ['lswi']);
  },
  
  applyMSAVI: function(image) {
    // apply MSAVI to an image
    var msavi = image.expression("(2*b('nir')+1-sqrt(pow((2*b('nir')+1),2)-8*(b('nir')-b('red'))))/2");
    return msavi.select([0], ['msavi']);
  },
  
  applySB: function(image) {
    // apply Soil Brightness to an image
    var sb = image.expression("0.3029*b('blue')+0.2786*b('green')+0.4733*b('red')+0.5599*b('nir')+0.5080**b('swir1')+0.1872*b('swir2')");
    return sb.select([0], ['sb']);
  },
  
  applyWRI: function(image) {
    // apply WRI to an image
    var wri = image.expression("(b('green')+b('red'))/(b('nir')+b('swir2'))");
    return wri.select([0], ['wri']);
  },
  
  applySWIR1: function(image) {
    // apply SWIR to an image
    var swir1 = image.expression("b('swir1')");
    return swir1.select([0], ['swir1']);
  },
  
  applyNIR: function(image) {
    // apply NIR to an image
    var nir = image.expression("b('nir')");
    return nir.select([0], ['nir']);
  },
  
  applyRED: function(image) {
    // apply RED to an image
    var red = image.expression("b('red')");
    return red.select([0], ['red']);
  },
  
    applyGREEN: function(image) {
    // apply GREEN to an image
    var green = image.expression("b('green')");
    return green.select([0], ['green']);
  },
  
    applyBLUE: function(image) {
    // apply BLUE to an image
    var blue = image.expression("b('blue')");
    return blue.select([0], ['blue']);
  },
  
    applySWIR2: function(image) {
    // apply BLUE to an image
    var swir2 = image.expression("b('swir2')");
    return swir2.select([0], ['swir2']);
  }
};

var getQABits = function(image, start, end, newName) {
    // Compute the bits we need to extract.
    var pattern = 0;
    for (var i = start; i <= end; i++) {
       pattern += Math.pow(2, i);
    }
    // Return a single band image of the extracted QA bits, giving the band
    // a new name.
    return image.select([0], [newName])
                  .bitwiseAnd(pattern)
                  .rightShift(start);
  };
  
  // A function to mask out cloudy shadow pixels.
var cloud_shadows = function(image) {
  // Select the QA band.
  var QA = image.select(['pixel_qa']);
  // Get the internal_cloud_algorithm_flag bit.
  return getQABits(QA, 3,3, 'Cloud_shadows').eq(0);
      // Return an image masking out cloudy areas.
  };
  
  // A function to mask out cloud pixels.
var clouds = function(image) {
  // Select the QA band.
  var QA = image.select(['pixel_qa']);
  // Get the internal_cloud_algorithm_flag bit.
  return getQABits(QA, 5,5, 'Cloud').eq(0);
  // Return an image masking out cloudy areas.
};

var maskClouds = function(image) {
  var cs = cloud_shadows(image);
  var c = clouds(image);
  image = image.updateMask(cs);
  return image.updateMask(c);
  };

// 2. Data Imports & Processing

// images
function generateLandsatCollection(){
  var L4collection = ee.ImageCollection('LANDSAT/LT04/C01/T1_SR')
      .filterDate(globOptions.startDate, globOptions.endDate)
      .map(maskClouds)
      .select(globOptions.bands7, globOptions.bandSelect);
  var L5collection = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
      .filterDate(globOptions.startDate, globOptions.endDate)
      .map(maskClouds)
      .select(globOptions.bands7, globOptions.bandSelect);
  var L7collection = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
      .filterDate(globOptions.startDate, globOptions.endDate)
      .map(maskClouds)
      .select(globOptions.bands7, globOptions.bandSelect);
  var L8collection = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
      .filterDate(globOptions.startDate, globOptions.endDate)
      .map(maskClouds)
      .select(globOptions.bands8, globOptions.bandSelect);
  var collectionFull = ee.ImageCollection(L4collection
      .merge(L5collection)
      .merge(L7collection)
      .merge(L8collection));
      //.filterBounds(site)
      //.filter(ee.Filter.intersects('.geo', globCoast.geometry(), null, null, 1000))
      //.filterMetadata('WRS_ROW', 'less_than', 120); 
  return collectionFull;
}
var collection = generateLandsatCollection();

//Add test
var non_mosaic = collection.filterBounds(spatialJoined_points);

//Will be used for frequency calculaton
var clouds_free = non_mosaic.map(maskClouds);

//display the mosaicked origional image and cloud free image
var mosaic_free = non_mosaic.map(maskClouds).median();
var visParams = {bands: ['red','green','blue'],min: [0,0,0],max: [2000, 2000, 2000]};
//Map.addLayer(non_mosaic, visParams, 'Merged Landsat Collection'); 

//Reducers

var reducer = ee.Reducer.min()
    .combine(ee.Reducer.max(), '', true)
    .combine(ee.Reducer.stdDev(), '', true)
    .combine(ee.Reducer.median(), '', true)
    .combine(ee.Reducer.intervalMean(25, 75).setOutputs(['int2575mean']), '', true)
    .combine(ee.Reducer.intervalMean(10, 90).setOutputs(['int1090mean']), '', true)
    //.combine(ee.Reducer.intervalMean(40, 60).setOutputs(['int4060mean']), '', true);
    //.combine(ee.Reducer.percentile([90]).setOutputs(['p75']), '', true);

// Data processing
var covariates = {
    aweiReduced: collection.map(landsatFunctions.applyAWEI)
        .reduce(reducer, globOptions.parallelScale),
    ndwiReduced: collection.map(landsatFunctions.applyNDWI)
        .reduce(reducer, globOptions.parallelScale),
    mndwiReduced: collection.map(landsatFunctions.applyMNDWI)
        .reduce(reducer, globOptions.parallelScale),
    ndviReduced: collection.map(landsatFunctions.applyNDVI)
       .reduce(reducer, globOptions.parallelScale),
    eviReduced: collection.map(landsatFunctions.applyEVI)
       .reduce(reducer, globOptions.parallelScale),
    lswiReduced: collection.map(landsatFunctions.applyLSWI)
       .reduce(reducer, globOptions.parallelScale),
    msaviReduced: collection.map(landsatFunctions.applyMSAVI)
       .reduce(reducer, globOptions.parallelScale),
    sbReduced: collection.map(landsatFunctions.applySB)
       .reduce(reducer, globOptions.parallelScale),
    wriReduced: collection.map(landsatFunctions.applyWRI)
       .reduce(reducer, globOptions.parallelScale),
    redReduced: collection.map(landsatFunctions.applyRED)
       .reduce(reducer, globOptions.parallelScale),
    greenReduced: collection.map(landsatFunctions.applyGREEN)
       .reduce(reducer, globOptions.parallelScale),
    blueReduced: collection.map(landsatFunctions.applyBLUE)
       .reduce(reducer, globOptions.parallelScale),
    swir1Reduced: collection.map(landsatFunctions.applySWIR1)
       .reduce(reducer, globOptions.parallelScale),
    swir2Reduced: collection.map(landsatFunctions.applySWIR2)
       .reduce(reducer, globOptions.parallelScale),
    nirReduced: collection.map(landsatFunctions.applyNIR)
       .reduce(reducer, globOptions.parallelScale),
    etopo: ee.Image('NOAA/NGDC/ETOPO1')
        .select(['bedrock'], ['etopo'])
        .resample('bicubic'),
    swOccurrence: ee.Image('JRC/GSW1_0/GlobalSurfaceWater')
        .select(['occurrence'], ['swOccurrence'])
        .unmask()
};
var trainComposite = covariates.aweiReduced
      .addBands(covariates.ndwiReduced)
//    .addBands(covariates.mndwiReduced)
//    .addBands(covariates.ndviReduced)
    .addBands(covariates.eviReduced)
//   .addBands(covariates.lswiReduced)
//    .addBands(covariates.msaviReduced)
    .addBands(covariates.sbReduced)
//    .addBands(covariates.wriReduced)
//    .addBands(covariates.redReduced)
//    .addBands(covariates.greenReduced)
//    .addBands(covariates.blueReduced)
//    .addBands(covariates.nirReduced)
//    .addBands(covariates.swir1Reduced)
     .addBands(covariates.swir2Reduced)
//    .addBands(covariates.etopo)
//    .addBands(covariates.swOccurrence);

var bands = trainComposite.bandNames();

// 3. Masking
var coastMask = US_coast
    .distance(globOptions.maskDistance).gte(-20); 
//Map.addLayer(coastMask, {color: 'FF00FF'}, 'coastMask');
//var coastMask2 = Regional_coast
//   .distance(globOptions.maskDistance).gte(2000); 
//Map.addLayer(coastMask2, {color: '00FF00'}, 'coastMask2');
var topoMask = covariates.etopo
    .gte(globOptions.maskDepth)
    .and(covariates.etopo.lte(globOptions.maskAltitude));
var topoMask = topoMask.updateMask(topoMask);
if (globOptions.maskApplySRTM) {
  var srtm = ee.Image('USGS/SRTMGL1_003')
    .lte(0);
  var finalMask = coastMask.multiply(topoMask)
    .rename('datamask')
    .byte()
    .updateMask(srtm);
} else {
  var finalMask = coastMask.multiply(topoMask)
    .rename('datamask')
    .byte();
}


// Overlay the points on the imagery to get training.
var full_points = trainComposite.select(bands).sampleRegions({
  collection: spatialJoined_points,
  properties: ['CLASS'],
  scale: 30
});

var withRandom = full_points.randomColumn('random');
//print('Full Point Set',withRandom);

var trainingPartition = withRandom.filter(ee.Filter.lt('random', globOptions.trainingValidationRatio));
var testingPartition = withRandom.filter(ee.Filter.gte('random', globOptions.trainingValidationRatio));

//print('Training Set:',trainingPartition);
//print('Validation Set:',testingPartition);

// 5. Classify 
var classifier = ee.Classifier.smileRandomForest({
    numberOfTrees: globOptions.nTrees, 
    //variablesPerSplit: 0,
    bagFraction: 0.5,
    seed: 0})
  .train(trainingPartition, 'CLASS', bands)
  .setOutputMode('CLASSIFICATION');
var classified = trainComposite.select(bands)
  .classify(classifier);
var finalOut = classified.byte()
  .mask(finalMask); 

// Postprocess
 var finalOut = finalOut.mask(finalOut
  .connectedPixelCount(globOptions.conPixels)
  .gte(globOptions.conPixels)); 
 
var TF_result = finalOut.updateMask(finalOut.eq(2)); // tf only


// Extra post-process

//var coastMask5k = ee.Image(1).clip(terrestrial5k);
//var invcoastMask5k = coastMask5k.mask(coastMask5k.mask().not());
//var finalOut = finalOut.updateMask(invcoastMask5k);

// 8. Run
var finalOutViz = {min: 0, max: 1, palette: ['0000FF', 'FF0000']};
Map.setCenter(-88.154331,30.672508, 7);
Map.addLayer(TF_result, finalOutViz, 'Tidal Flats',false);
Export.image.toDrive({
  image: TF_result.clip(aoi),
  description: '2020_tidal_flats_MS_AL_LA_new0215',
  scale: 30,
  region: TF_result.clip(aoi),
  maxPixels: 10000000000000
});


//Compare the landcover of your validation data against the classification result
//var testAccuracy = validation.errorMatrix('CLASS', 'classification');
//Print the error matrix to the console
//print('Validation error matrix: ', testAccuracy);
//Print the overall accuracy to the console
//print('Validation overall accuracy: ', testAccuracy.accuracy());

// Classify the test FeatureCollection.
var test = testingPartition.classify(classifier);

// Print the confusion matrix.
var confusionMatrix = test.errorMatrix('CLASS', 'classification');
//print('Confusion Matrix', confusionMatrix);
print('Validation overall accuracy: ', confusionMatrix.accuracy());
