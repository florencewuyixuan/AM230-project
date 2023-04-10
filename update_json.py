import json

def update_json(filename, velocity):
    jsonData = json.load(open(filename, "r"))
    num = len(velocity)
    for i in range(num):
        jsonData['system']['particles'][i]["v"] = velocity[i]*jsonData['system']['particles'][i]["n"]

    json.dump(jsonData, open(filename ,"w"), indent=4)
