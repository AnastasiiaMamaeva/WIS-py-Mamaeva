import numpy as np
import random


    

correct = False
again = True
while again == True: 
    random_number = random.randint(1, 20)
    i = 1
    while correct == False:
      user_num = int(input("try to guess a number between 1 and 20"))
      if random_number < user_num:
          print("Oh, your number is too big, try again!")
      elif random_number > user_num:
          print("Oh, your number is too small, try again!")
      elif random_number == user_num:
          print("Good job! You guessed right! You did", i, "tries!", sep = " ")
          correct = True
          ask = input("Do you want to play again? (type 'Yes' or 'No')")
          if ask == "Yes":
              break
          else:
              again == False
      i +=1