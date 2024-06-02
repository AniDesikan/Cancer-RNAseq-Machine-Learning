#!/usr/bin/python

import hashlib
import pymysql
import sys
import cgi
import cgitb
import json
import csv
import os
import io

cgitb.enable()
####################
## FUNCTIONS ##
####################

######
# Function to connect to the localhost database
def connect_db():
        #connect to database
        connection = pymysql.connect(host="localhost", user="root", password="6HcCzSCZVhXrfvK!", db="hotresdb", port=3306)
        #get cursor
        cursor = connection.cursor()
        return cursor, connection

#######
# Function to disconnect from the localhost database
def disconnect_db(cursor, connection):
        #close cursor and connection
        cursor.close()
        connection.close()


#######
# Function to get the EID from a study name in the EXPERIMENT table

def get_EID(study):
        query = """
        SELECT EID
        FROM EXPERIMENT_TEST
        WHERE STUDY_NAME = %s
        """
        try:
                cursor.execute(query, (study,))
                result = cursor.fetchone()
                if result:
                        return result[0]
                else:
                        return None
        except pymysql.Error as e: 
                print(e, query)
                return None

#######
# Function to parse the samples csv to input the results into the samples table

def parse_sample_csv(csv_path):
        with open(csv_path, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                        # Extract data from each row
                        sample_id = row['sample_id']
                        DPI = row['DPI']
                        virus_count = row['virus_count']
                        Experimental_condition = row['Experimental_condition']
                        PID = row['PID']
                return sample_id, DPI, virus_count, Experimental_condition, PID

###########################################

print("Content-type: text/html\n")

content_length = int(os.environ.get('CONTENT_LENGTH', 0))
raw_data = sys.stdin.buffer.read(content_length)

# Parse the raw request body as FormData
form = cgi.FieldStorage(fp=(raw_data), environ=os.environ, keep_blank_values=True)

# form = cgi.FieldStorage()
if (form):
        cursor, connection = connect_db()

        # This selector is what chooses which query is executed, chosen in the HTML file by which button is pressed
        selector = form.getvalue("selector", "")
        
        if (selector == "papers"):
                query = """
                select PUBID
                from PUBLICATION
                order by PUBID asc
                """
                try:
                        cursor.execute(query)	
                except pymysql.Error as e: 
                        print(e,query)	
        
        elif (selector == "experiment"):
                query = """
                select STUDY_NAME
                from EXPERIMENT_TEST
                order by STUDY_NAME asc
                """
                try:
                        cursor.execute(query)	
                except pymysql.Error as e: 
                        print(e,query)	
        
        elif (selector == "enter_experiment"):
                study_name = form.getvalue("studyname")
                virus = form.getvalue("virus")
                strain = form.getvalue("strain")
                species = form.getvalue("species")
                tissue_type = form.getvalue("tissuetype")
                exposure_route = form.getvalue("exposureroute")
                pubid = form.getvalue("paper")

                query = """
                INSERT INTO EXPERIMENT_TEST (STUDY_NAME, VIRUS, STRAIN, SPECIES, TISSUE_TYPE, EXPOSURE_ROUTE, SEQUENCING_TYPE, PUBID)
                VALUES (%s, %s, %s, %s, %s, %s, "RNAseq", %s)
                """

                try:
                        cursor.execute(query, (study_name, virus, strain, species, tissue_type, exposure_route, pubid))
                        connection.commit()
                except pymysql.Error as e:
                        response = {"success": False, "message": str(e)}
                finally:
                        disconnect_db(cursor, connection)

        elif (selector == "enter_samples"):

                study_name = form.getvalue("experiment")
                sample_csv = form.getvalue("samplecsv")
                EID = get_EID(study_name)
                sample_id, DPI, virus_count, Experimental_condition, PID = parse_sample_csv(sample_csv)

                query = """
                INSERT INTO SAMPLES_TEST (sample_id, DPI, virus_count, Experimental_condition, PID, EID)
                VALUES (%s, %s, %s, %s, %s, %s)
                """
                try:
                        cursor.execute(query, (sample_id, DPI, virus_count, Experimental_condition, PID, EID))
                        connection.commit()
                except pymysql.Error as e:
                        response = {"success": False, "message": str(e)}
                finally:
                        disconnect_db(cursor, connection)

        results = cursor.fetchall()


print(json.dumps(results))