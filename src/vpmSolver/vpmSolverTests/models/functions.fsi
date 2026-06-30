&TRIAD
  id = 25
  extId = 1
  extDescr = 'T1'
  nDOFs = 6
  ur = 1.0 0.0 0.0 0.0
       0.0 1.0 0.0 0.0
       0.0 0.0 1.0 0.0
/

&SENSOR
  id = 16
  extId = 1
  extDescr = 'Time'
  type = 'TIME'
/

&SENSOR
  id = 2611
  extId = 1
  extDescr = 'Sensor on T1.x'
  type = 'TRIAD'
  triadId = 25
  dof = 1
  dofEntity = 'POS'
  dofSystem = 'GLOBAL'
/

&SENSOR
  id = 2621
  extId = 1
  extDescr = 'Sensor on T1.y'
  type = 'TRIAD'
  triadId = 25
  dof = 2
  dofEntity = 'POS'
  dofSystem = 'GLOBAL'
/

&ENGINE
  id = 17
  extId = 1
  extDescr = 'f1'
  functionId = 18
  nArg = 1, argSensorId = 16
/

&ENGINE
  id = 19
  extId = 2
  extDescr = 'f2'
  functionId = 20
  nArg = 1, argSensorId = 16
/

&ENGINE
  id = 21
  extId = 3
  extDescr = 'f3(x,y,z)'
  functionId = 22
  nArg = 3, argSensorId = 2611 2621 16
/

&ENGINE
  id = 23
  extId = 4
  extDescr = 'f4'
  functionId = 24
  nArg = 1, argSensorId = 16
/

&FUNCTION
  id = 18
  extId = 1
  type = 'SINUSOIDAL'
  realDataSize = 5
  realData = 1.0 0.0 1.0 0.0 0.0
/

&FUNCTION
  id = 24
  extId = 1
  type = 'SCALE'
  realDataSize = 1
  realData = 1.0
/

&FUNCTION
  id = 22
  extId = 1
  type = 'MATH_EXPRESSION'
  expression = 'x+y*y-xyz'
  nArg = 3
/

&FUNCTION
  id = 20
  extId = 1
  type = 'RAMP'
  realDataSize = 3
  realData = 0.0 1.0 0.0
/
