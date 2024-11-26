import numpy as np
import random


    

def check_input(user_input, random_number, attempts):
    if random_number < int(user_input):
        print("Oh, your number is too big, try again!")
    elif random_number > int(user_input):
        print("Oh, your number is too small, try again!")
    elif random_number == int(user_input):
        print("Good job! You guessed right! You did", attempts, "tries!", sep=" ")
        return True
    return False


def main():
    correct = False
    random_number = random.randint(1, 20)
    attempts = 1

    while not correct:
        user_input = input("Try to guess a number between 1 and 20: ")
        if user_input == "x":
            return False
        elif user_input == "s":
            print("The number is", random_number, sep=" ")
        elif user_input == "n":
            if input("Do you want to play a new game? (Yes/No): ").lower() == "yes":
                return True
        elif user_input.isdigit():
            correct = check_input(user_input, random_number, attempts)
            if correct:
                if input("Do you want to play again? (type 'Yes' or 'No'): ").lower() != "yes":
                    return False
        else:
            print("Invalid input. Please enter a number, 's', 'n', or 'x'.")
        attempts += 1


def play_game():
    play_again = True
    while play_again:
        play_again = main()


if __name__ == "__main__":
    play_game()
