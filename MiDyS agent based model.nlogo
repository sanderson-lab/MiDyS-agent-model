globals
[
  tick-counter     ;; counts number of ticks
  FusC             ;; fusion counter
  FisC             ;; fission counter
  FFI              ;; fusion/fission index
  locals           ;; temp agentset candidates to fuse with
]

turtles-own
[
  age              ;; number of ticks turtle has existed
  fusion           ;; probability of fusion (average)
  fusion-pot       ;; probability of fusion with membrane potential variable
  fusion-mt        ;; probability of fusion with mitotimer variable
  fission          ;; probability of fission (average)
  fission-pot      ;; probability of fission with membrane potential variable
  fission-mt       ;; probability of fission with mitotimer variable
  links-count      ;; number of links for each turtle
  area             ;; area of mitochondrial object
  potential        ;; mitochondrial membrance potential
  mitotimer        ;; mitotimer ratio
  major            ;; longest axis of mitochondrial object
  minor            ;; shortest axis of mitochondrial object
  AR               ;; major / minor
  object-num       ;; label for mitochondrial object
  numbered?        ;; check during renumbering
]

links-own
[
  link-num
]

;;;;;;;;;;;;;;;;;;;;;;;;
;;; Setup Procedures ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all                                            ;; clear previous iterations of model
  set-default-shape turtles "circle"                   ;; create first configuration of turtles
  create-turtles Initial_Network_Density               ;; create turtles with random headings
  set tick-counter 1                                   ;; start tick timer
  set FusC 0                                           ;; intial setting of FusC
  set FisC 0                                           ;; intial setting of FusC
  set FFI 0                                            ;; initial setting of FFI
  ask turtles
  [
    set color blue                                     ;; set starting parameters for all turtles
    set size .5
    set potential random-normal -0.03418 0.2121        ;; assign membrane to gaussian distribution
    set age 0
    set mitotimer random-normal -0.5969 0.1825         ;; assign mitotimer to gaussian distribution
    set xcor random-xcor                               ;; move turtles to random starting location
    set ycor random-ycor
    while [any? other turtles in-radius 1][rt random-float 360 forward random-float 1]   ; move turtles to remove overlap
    set object-num who                                 ;; assign initial object numbers
  ]
  ask turtles
  [
    let z count other turtles in-radius 2 with [not link-neighbor? myself]                                 ;; agents make starting networks
    while [z >= 1] [
      if random 100 <= initial-networking [                                                                  ;; initial networking variable determines fusion probability
        create-link-with one-of other turtles in-radius 2 with [not link-neighbor? myself]
        if object-num > (min [object-num] of link-neighbors) [set object-num (min [object-num] of link-neighbors)]
        ask link-neighbors [set object-num [object-num] of myself]
      ]
      set z z / 2
    ]
    set links-count count link-neighbors

    if links-count > 0
    [
      ifelse links-count >= max [links-count] of link-neighbors                 ;; set object number to be equal across links
    [
      ask link-neighbors
      [
        set object-num [object-num] of myself
        ask my-links
        [
         set link-num [object-num] of myself
        set thickness 0.3
        ]
      ]
    ]
      [
        set object-num min [object-num] of link-neighbors
        ask my-links
        [
         set link-num [object-num] of myself
         set thickness 0.3
        ]
      ]
   ]
  ]
  renumber                                  ;; renumbers objects
  manipulation                              ;; execute membrane potential or oxidation manipulation
  ask turtles [calc]
  ask turtles [calc2]
  update-view
  reset-ticks
end

;;;;;;;;;;;;;;;;;;;;;;;
;;; Main Procedures ;;;
;;;;;;;;;;;;;;;;;;;;;;;

to go
  set FusC 0                                           ;; reset of FusC
  set FisC 0                                           ;; reset of FusC
  manipulation                                             ;; execute membrane potential or oxidation manipulation
  ask turtles                                                                              ;; behavior command for turtles
  [
    migrate
    calc
  ]
  ask turtles
  [
    calc2
    fiss
  ]
  renumber
  ask turtles
  [
    calc
  ]
  ask turtles
  [
    calc2
    fuse
  ]
  renumber
  ask turtles
  [
    mitophagy
    aging
  ]
  repeat (random-float 0.05 * (count turtles)) [biogenesis]        ;; growth of new agents
  update-view
  set FFI (((count turtles with [count link-neighbors > 1]) / (count turtles)) * 2) - 1
  set tick-counter tick-counter + 1
  tick
end


;; functions ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to migrate                                                                                                             ;; movement function
  if random 100 > 50
  [
  ifelse count my-links > 0
  [
    if (max [link-length] of my-links <= 1)               ;; agents cannot move too far from other connected agents
    [
      rt random 360
      forward random-float 0.5
      while [any? other turtles in-radius 1][rt random 360 forward random-float 0.5]        ;; remove overlap
    ]
  ]
  [
    rt random 360
    forward random-float 0.5
    while [any? other turtles in-radius 1][rt random 360 forward random-float 0.5]
  ]
  ]
  if count my-links > 0 [
  if [link-length] of max-one-of my-links [link-length] > 4 [
    ask my-links with [link-length > 4]
    [die]
  ]]
end

to calc                                                                   ;; calculate variables for each individual agent
    set links-count count link-neighbors                                  ;; count link neighbors
    set area ((count (turtles with [object-num = ([object-num] of myself)])) * 0.5) + (0.3 * sum [link-length] of links with [link-num = [object-num] of myself])        ;; area = (0.5 * connected agents and self) + (0.5 * sum of link length)
    ifelse links-count > 0 [set major (distance max-one-of other turtles with [object-num = ([object-num] of myself)][distance myself])][set major 0.5]                  ;; major = distance to farthest link neighbor
    if major < 0.5 [set major 0.5]
    ifelse links-count > 0 [set minor (distance min-one-of other turtles with [object-num = ([object-num] of myself)][distance myself])][set minor 0.3]                  ;; minor = distance to nearest link neighbor
    if minor < 0.3 [set minor 0.5]
end

to calc2                  ;; equalize variables across agents in same object
    set major max [major] of turtles with [object-num = ([object-num] of myself)]     ;; set major
    set minor min [minor] of turtles with [object-num = ([object-num] of myself)]     ;; set minor
    set AR (major / minor)                                                            ;; set aspect ratio
    set fusion-pot 100 * ((e ^ (-2.452 + ((log area 10) * -3.410) + (major * 0.1747) + (minor * 0.8830) + ((log AR 10) * 1.939) + (potential * 0.7833))) / (1 + (e ^ (-2.452 + ((log area 10) * -3.410) + (major * 0.1747) + (minor * 0.8830) + ((log AR 10) * 1.939) + (potential * 0.7833)))))
    set fusion-mt 100 * ((e ^ (-5.447 + ((log area 10) * -5.938) + (major * 0.5648) + (minor * 4.738) + ((log AR 10) * 2.413) + (mitotimer * 2.599))) / (1 + (e ^ (-5.447 + ((log area 10) * -5.938) + (major * 0.5648) + (minor * 4.738) + ((log AR 10) * 2.413) + (mitotimer * 2.599)))))
    if fusion-pot > fusion-mt [set fusion fusion-pot]
    if fusion-mt > fusion-pot [set fusion fusion-mt]
    set fission-pot 100 * ((e ^ (-2.688 + ((log area 10) * 3.852) + (major * -0.2560) + (minor * -0.06415) + ((log AR 10) * 2.027) + (potential * -1.485))) / (1 + (e ^ (-2.688 + ((log area 10) * 3.852) + (major * -0.2560) + (minor * -0.06415) + ((log AR 10) * 2.027) + (potential * -1.485)))))
    set fission-mt 100 * ((e ^ (-6.048 + ((log area 10) * 3.712) + (major * -0.3736) + (minor * 1.531) + ((log AR 10) * 4.284) + (mitotimer * -0.9781))) / (1 + (e ^ (-6.048 + ((log area 10) * 3.712) + (major * -0.3736) + (minor * 1.531) + ((log AR 10) * 4.284) + (mitotimer * -0.9781)))))
    if fission-pot > fission-mt [set fission fission-pot]
    if fission-mt > fission-pot [set fission fission-mt]
    set FFI (((count turtles with [count link-neighbors > 1]) / (count turtles)) * 2) - 1
end


to fuse                    ;; fusion function
  set locals other turtles in-radius 2.5                        ;; find fusion candidates
  set locals locals with [object-num != [object-num] of myself] ;; remove fusion candidates already in object
  if count locals > 0 [
    if random 100 < fusion
    [
      create-link-with max-one-of locals [fusion] [set thickness 0.3]      ;; fuse with one fusion candidate
      set object-num min [object-num] of link-neighbors
      ask my-links [set link-num [object-num] of myself]
      set FusC FusC + 1
    ]
   ]
end

to renumber
  ask turtles [
    set numbered? false
    set links-count count link-neighbors
  ]
  loop
  [
    let start min-one-of turtles with [not numbered?] [object-num]
    if start = nobody [ stop ]
    ask start [
    if links-count > 0
    [
        connect (min [object-num] of link-neighbors)
    ]
    if links-count = 0
    [
      set numbered? true
      if count other turtles with [object-num = [object-num] of myself] > 0
      [
        set object-num (max [object-num] of turtles + 1)
        ]
    ]
    ]
  ]
end

to connect [object]
  if numbered? [ stop ]
  set numbered? true
  set object-num object
  ask my-links [set link-num object]
  ask link-neighbors [connect object]
end


to fiss                             ;; fission function
  if count my-links > 0
  [
    if random 100 < fission
    [ask link-with max-one-of link-neighbors [fission] [die]      ;; fission of one link
      set object-num ((max [object-num] of turtles) + 1)
      set FisC FisC + 1
    ]
  ]

end

to aging                                                                                ;; age turtles and alter variables
  set age age + 1
  set potential potential + (random-float 0.07)                                         ;; membrane potential flucuates in both directions
  set potential potential - (random-float 0.07)
  set mitotimer mitotimer - (random-float 0.05 * (abs potential))                         ;; mitotimer shifts negative over time with influence by membrane potential
end

to biogenesis
  if count turtles with [potential >= -0.2 AND potential <= 0.2 AND count turtles-on neighbors < 4] > 0
  [
  ask one-of turtles with [potential >= -0.2 AND potential <= 0.2 AND count turtles-on neighbors < 4]              ;; selects healthy turtle to copy
  [hatch 1                                                                                                         ;; create copy and set parameters
    [
      set age 0
      set potential random-normal -0.03418 0.2121
      set mitotimer mitotimer + 0.14
      migrate
    ]
    create-links-with other turtles with [age = 0] in-radius 2 [set thickness 0.3]                               ;; original turtle fuses with copy
  ]
  ]
end

to mitophagy                                                 ;; mitophagy to remove turtles
   if links-count < 1
   [
     if random 100 > 90                                      ;; 10% chance that single turtle will be removed
      [ask my-links [die]
        die]
   ]
end

to update-view
  if view-option = "Default"
  [
    ask turtles [set color 95]
    ask links [set color 96]
  ]
  if view-option = "Membrane Potential"
 [
    ask turtles [set color (potential * 5) + 15]
    ask links
    [
    set color ([color] of end1 + [color] of end2) / 2
    ]
  ]
  if view-option = "MitoTimer"
  [
    ask turtles [set color (mitotimer * 6.42) + 89.3]
    ask links
    [
      set color ([color] of end1 + [color] of end2) / 2
    ]
  ]
  if view-option = "Objects"
  [
     ask turtles [set color (object-num + 0.2)]
     ask links
    [
      set color ([color] of end1)
    ]
  ]
end

to manipulation
    ask turtles [
    set potential potential + (0.1 * (Polarization / 100))
    set mitotimer mitotimer - (0.1 * (Oxidation / 100))
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
204
10
919
726
-1
-1
7.0
1
10
1
1
1
0
0
0
1
-50
50
-50
50
1
1
1
ticks
30.0

BUTTON
15
10
197
45
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
17
50
198
83
Initial_Network_Density
Initial_Network_Density
100
2000
1000.0
100
1
NIL
HORIZONTAL

SLIDER
20
89
190
122
initial-networking
initial-networking
0
100
70.0
1
1
NIL
HORIZONTAL

BUTTON
15
186
195
219
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
15
224
195
259
go once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
925
10
1125
160
Age Distribution
Age
Count
0.0
10.0
0.0
10.0
true
false
"set-plot-x-range 0 10\nset-plot-y-range 0 count turtles\nset-histogram-num-bars 7" "set-plot-x-range 0 (max [age] of turtles + 1)\nset-plot-y-range 0 count turtles\nset-histogram-num-bars 7"
PENS
"Age" 1.0 1 -16777216 true "" "histogram [age] of turtles"

PLOT
926
164
1126
314
# of Mitochondria
Time
Mitochondria
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"Mitochondria" 1.0 0 -16777216 true "" "plot count turtles"

PLOT
1130
475
1330
624
Fission/Fusion
Time
Events
0.0
10.0
0.0
100.0
true
false
"" ""
PENS
"FusC" 1.0 0 -16777216 true "" "plot (FusC / count turtles * 100)"
"FisC" 1.0 0 -2674135 true "" "plot (FisC / count turtles * 100)"

PLOT
1129
10
1329
160
Average Membrane Potential
Potential
Count
0.0
10.0
0.0
10.0
true
false
"set-plot-y-range -1 1" ""
PENS
"potential" 1.0 0 -16777216 true "" "plot mean [potential] of turtles"

PLOT
926
319
1126
469
Fusion Fission Probabilities
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"p-Fusion" 1.0 0 -16777216 true "" "plot mean [fusion] of turtles"
"p-Fission" 1.0 0 -5298144 true "" "plot mean [fission] of turtles"

CHOOSER
25
131
183
176
View-Option
View-Option
"Default" "Membrane Potential" "MitoTimer" "Objects"
0

PLOT
1334
10
1534
160
Membrane Potential
Membrane Potential
Count
0.0
10.0
0.0
10.0
true
false
"set-plot-x-range -1 1\nset-histogram-num-bars 20" ""
PENS
"potential" 1.0 0 -16777216 true "" "histogram [potential] of turtles"

PLOT
1336
167
1536
317
MitoTimer
MitoTimer
Count
0.0
10.0
0.0
10.0
true
false
"set-plot-x-range -1.2 .2\nset-histogram-num-bars 20" ""
PENS
"mitotimer" 1.0 0 -16777216 true "" "histogram [mitotimer] of turtles"

PLOT
1130
320
1330
470
Area Distribution
Area
Count
0.0
10.0
0.0
10.0
true
false
"set-plot-x-range 0 max [area] of turtles\nset-histogram-num-bars 20" ""
PENS
"area" 1.0 1 -16777216 true "" "histogram [area] of turtles"

PLOT
927
475
1127
625
Fusion Probability
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"set-plot-x-range 0 100\nset-histogram-num-bars 20" ""
PENS
"p-fusion" 1.0 0 -16777216 true "" "histogram [fusion] of turtles"

PLOT
927
629
1127
779
Fission Probability
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"set-plot-x-range 0 100\nset-histogram-num-bars 20" ""
PENS
"p-fission" 1.0 0 -16777216 true "" "histogram [fission] of turtles"

PLOT
1131
166
1331
316
Average MitoTimer
MitoTimer
Count
0.0
10.0
0.0
10.0
true
false
"set-plot-y-range -1.2 0.2\n" ""
PENS
"mitotimer" 1.0 0 -16777216 true "" "plot mean [mitotimer] of turtles"

SLIDER
19
266
191
299
Polarization
Polarization
-100
100
0.0
10
1
%
HORIZONTAL

PLOT
1132
630
1329
780
FFI
Time
FFI
0.0
10.0
0.0
10.0
true
false
"set-plot-y-range -1 1" ""
PENS
"ffi" 1.0 0 -16777216 true "" "plot FFI"

SLIDER
19
305
191
338
Oxidation
Oxidation
0
100
0.0
10
1
%
HORIZONTAL

BUTTON
18
345
190
378
Export All Plot Data
export-all-plots Export_File
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
18
386
191
446
Export_File
plots.csv
1
0
String

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
