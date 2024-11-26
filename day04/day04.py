import numpy as np
import random


    

correct = False
again = True
while again == True: 
    random_number = random.randint(1, 20)
    i = 1
    while correct == False:
      user_num = input("try to guess a number between 1 and 20")
      if user_num == "x":
        again = False
        break
      elif user_num == "n":
        if input("Do you want to play a new gaim?") == "Yes":
          break
      elif user_num == "s":
          print("The numbr is", random_number, sep = " ")
      elif random_number < int(user_num):
          print("Oh, your number is too big, try again!")
      elif random_number > int(user_num):
          print("Oh, your number is too small, try again!")
      elif random_number == int(user_num):
          print("Good job! You guessed right! You did", i, "tries!", sep = " ")
          correct = True
          ask = input("Do you want to play again? (type 'Yes' or 'No')")
          if ask == "Yes":
            break
          elif ask == "x":
              again = False
              break
          elif ask == "n":
            if input("Do you want to play a new gaim?") == "Yes":
              break
          elif ask == "s":
            print("The numbr is", random_number, sep = " ")
          else:
            again = False
      i +=1