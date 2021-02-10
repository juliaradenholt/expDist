import expected_distance as e
import models as m
import random

models = { 'Dayhoff':m.Dayhoff(), 'JTT':m.JTT(), 'LG':m.LG(), 'WAG':m.WAG() }
case= ['W', 'X', 'Y', 'Z']
models_str = ['PAM', 'JTT', 'LG', 'WAG']
test_case = (random.randint(0, 3))
test_num = (random.randint(1, 35))
test_mod = (random.randint(0, 3))

print("\n \nThis program allows you to enter a substiution model and a filepath.")
print("The program will then estimate and print the expected values of the evolutionary distances")
print("between the sequences in the file. The only valid format for the files is the FASTA format.\n \n")

print("Enter path to file, if nothing is entered, a random file from the testcases is picked:")
filename = input()
print("Which substiution model do you want to use? Please enter: Dayhoff, JTT, LG or WAG:")
model_name = input()

try:
    model = (models.get(model_name))
    filepath = filename if filename != '' else ("./testcases/"+case[test_case]+"/"+str(test_num)+'/m'+models_str[test_case])+".fasta"
    filepath = filename if filename != '' else ("./testcases/"+case[test_case]+"/"+str(test_num)+'/m'+models_str[test_case])+".fasta"
    print("Running "+str(filepath)+" with "+str(model_name))
    expected_value = e.expected_dist(filepath, model)
except AttributeError or TypeError:
    print("Please enter a valid model and filepath.")
