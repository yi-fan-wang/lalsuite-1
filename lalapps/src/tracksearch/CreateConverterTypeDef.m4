define(`SCALE',0)
ifelse(TYPECODE,`D', `define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S', `define(`TYPE',`REAL4')')
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')
ifelse(TYPECODE,`',  `define(`TYPE',`REAL4')')
ifelse(TYPECODE,`D', `define(`SCALE',20)')

define(`XFUNC',`format(`XLALFrGet%sFrameConvertToREAL4TimeSeries',TYPE)')
define(`CREATESERIES',`format(`XLALCreate%sTimeSeries',TYPE)')
define(`GETMETA',`format(`XLALFrGet%sTimeSeriesMetaData',TYPE)')
define('GETDATA',`format(`XLALFrGet%sTimeSeriesData',TYPE)')
define(`DESTROYSERIES',`format(`XLALDestroy%sTimeSeries',TYPE)')
define(`VARTYPE',`format(`%sTimeSeries',TYPE)')

int XFUNC (REAL4TimeSeries *inputSeries FrStream *stream);
