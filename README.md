Pipeline's desined in Model–View–Controller architecture. Main parts of the script based on Match script, which you should have installed in order to have working pipeline. More information you can find here: http://spiff.rit.edu/match/match-0.16/match.html

1) For start create folder named with name of your cluster exp. "Gulliver_27/". 

2) In that folder create .yaml file with name of your cluster exp. "Gulliver_27.yaml"

3) In the same folder create new folder "input/" and put there your UV and GAIA data exp. "uv.csv" and "star_table_Gulliver_27.csv". 

4) To obtain GAIA data use "get_stars_at_oc.py". To use it type in command line:

   "python3 get_stars_at_oc.py Gulliver_27 146.0801 -54.1154 0.12 19"
   
   where: Gulliver_27 is the name of the object (We need it only to navigate in folders, in general, GAIA doesn't provied cluster names. 

   Be free to name it howewer you want, just make it the same name as folder and .yaml file from step 1) and 2)), 146.0801 is the RA, -54.1154 is the DEC, 0.12 is the search radius in degrees (Pick radius not too large, so it will be approxomatley the same as UV data       radius. Or you have a risk to have 0 matched stars.), and 19 is the magnitude limit. 

5) UV data should contain 3 columns with names: "ra", "dec", "u". First two stands for coordinates and last stand for magnitude in U-filter

   GAIA data obtained by "get_stars_at_oc.py" should be in the same form always. (Sometimes first colums SOURCE_ID could be in lowercase. If it's so, change it to uppercase)

6) In .yaml file add information about:
   your cluster: "cluster",

   coordinates of the center of GAIA data: "coords",

   magnitude you will use for matching: "gaia_mag",

   would you like to make plots: "plots" (True - make plots / False - don't make plots),

   parameters for matching: "matchrad", "trirad", "nobj",

   colours you want to use for CC-diagram: "colour1", "colour2"

   do you want to adjust positions of obtained stars on CC-diagram, in order to get they offset from MS in that filter: "adjust" (True - adjust / False - don't adjust)

7) Script has 4 procedures:
   
   "coords_reproject" takes RA and DEC coordinates of objects from UV and GAIA data and using project_coords from Match script to project them on tangential plane and get X and Y coordinates. More information avaliable by link above. 

   "matching" matching UV and GAIA stars and gives as output lists of matched stars from both datasets. 

   "filtering" filters stars from Hunt catalogue of clusters memberships with probability of membership >0.7. At the end you will have lists of filtered stars from both datasets. 

   "cc_diagram" creates CC-diagram in colours given by "colour1", "colour2" in .yaml file. If you have "adjust : True", firstly you will see interactive plot, when you can adjust vertical and horisontal positions of filtered stars (pink triangles). After you adjusted 
   them, so they will fit some part of MS, close poped up window and your plot will be saved in "plots/" directory. 

8) To run the script procedure by procedure, type in the command line:

   "python3 controller.py Gulliver_27 coords_reproject",
   
   where "Gulliver_27" is the name or your cluster and "coords_reproject" is the name of procedure.

11) In general case, you can run all the procedurer at ones by typing:

    "python3 controller.py Gulliver_27 process",

    where process is the name of procedure which contains all previously mentioned procedures. RECCOMENDED to do it in this way.  

14) All of the parameters should be changes in .yaml file. Try to avoid modifications of controller. and process_data.py if you're not sure what you're doing. 
