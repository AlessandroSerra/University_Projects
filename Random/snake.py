import pygame as pygame
import random 
import time

pygame.init()

dis_width = 600
dis_height = 800

dis = pygame.display.set_mode((dis_height, dis_width))

pygame.display.update()
pygame.display.set_caption('Snake game by AleSerra')

white = (255, 255, 255)
black = (0, 0, 0)
blue = (0, 0, 255)
red = (255, 0, 0)

snake_block = 10
snake_speed = 30

font_style = pygame.font.SysFont(None, 50)

clock = pygame.time.Clock()

def message(msg, color):

    mesg = font_style.render(msg, True, color)
    dis.blit(mesg, [dis_width / 2 , dis_width / 2])

def gameLoop():
    game_over = False
    game_close = False

    x1 = dis_width / 2
    y1 = dis_width / 2

    x1_change = 0
    y1_change = 0

    foodx = round(random.randrange(0, dis_width - snake_block) / 10.0) * 10.0
    foody = round(random.randrange(0, dis_width - snake_block) / 10.0) * 10.0

    while not game_over:

        while game_close == True:
            dis.fill(white)
            message("You Lost! Press Q-Quit or C-Play Again", red)
            pygame.display.update()

            for event in pygame.event.get():
                if event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_q:
                        game_over = True
                        game_close = False
                    if event.key == pygame.K_c:
                        gameLoop()

        for event in pygame.event.get():

            if event.type == pygame.QUIT:
                game_over = True

            if event.type == pygame.KEYDOWN:

                if event.key == pygame.K_LEFT:
                    x1_change = -10
                    y1_change = 0

                elif event.key == pygame.K_RIGHT:
                    x1_change = 10
                    y1_change = 0

                elif event.key == pygame.K_UP:
                    x1_change = 0
                    y1_change = -10

                elif event.key == pygame.K_DOWN:
                    x1_change = 0
                    y1_change = 10

        if x1 >= dis_width or x1 < 0 or y1 >= dis_height or y1 < 0:
            game_over = True

        x1 += x1_change
        y1 += y1_change

        dis.fill(white)
        pygame.draw.rect(dis, black, [x1, y1, 10, 10])
        pygame.draw.rect(dis, black, [x1, y1, snake_block, snake_block])

        pygame.display.update()

        if x1 == foodx and y1 == foody:
            print("Yummy!!")

        clock.tick(snake_speed)

    pygame.quit()

gameLoop()