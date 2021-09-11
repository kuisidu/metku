def log_to_file(function):
    def wrapper(*args, **kwargs):
        import sys
        import datetime as dt
        import inspect
        import os
        original_stdout = sys.stdout
        filename = f'{function.__name__}_{dt.datetime.now().strftime("%Y_%m_%d__%H_%M_%S")}'
        os.mkdir(f'{os.getcwd()}/{filename}')
        os.chdir(f'{os.getcwd()}/{filename}')
        print(f'writing to file: {filename}')
        with open(f'{filename}.txt', 'w', encoding='UTF-8') as file:
            sys.stdout = file

            print(inspect.getsource(sys.modules['__main__']))
            print('\n\n####################################### SOLUTION ##########################################\n\n')
            function(*args, **kwargs)

            sys.stdout = original_stdout
        os.chdir(f'{os.getcwd()}/..')
        print(f'done writing')
    return wrapper
